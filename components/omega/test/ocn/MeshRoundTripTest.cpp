//===-- Test driver for OMEGA mesh round trip ------------------*- C++ -*-===//
//
/// \file
/// \brief Test that the mesh read by Omega is written back out unchanged
///
/// This driver reads the horizontal mesh the way Omega does (Decomp and
/// HorzMesh, through the HorzMeshIn stream) and then writes the contents of
/// that stream back out through an Omega output stream. Every variable in the
/// output file is then compared with the variable of the same name (or its
/// older MPAS name) in the mesh file and must be identical, so that a stream
/// asked to write the horizontal mesh writes a valid MPAS mesh. Any
/// conversion Omega applies on reading must therefore be undone on writing.
/// Connectivity entries beyond the count of valid entries in each row are
/// padding and are not compared. The output must also carry the global
/// attributes HorzMesh requires of any mesh file (on_a_sphere, is_periodic
/// and, as applicable, sphere_radius or x_period and y_period) with the
/// input values.
///
/// Usage: testMeshRoundTrip.exe [MeshFile]
///
/// MeshFile is read by the HorzMeshIn stream (default OmegaMesh.nc).
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "mpi.h"
#include "pio.h"
#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace OMEGA;

namespace {

// Name of the output stream added to the configuration by this test
const std::string HorzOutStreamName = "MeshRoundTripHorzOut";

// Maximum number of mismatched entries to report for each variable
constexpr int MaxReport = 5;

//------------------------------------------------------------------------------
// Returns a write stream (OnDemand, so written only when forced) to OutFile
// with the given contents, based on a copy of an existing write stream
YAML::Node makeOutStream(const YAML::Node &Template, const std::string &OutFile,
                         const std::string &Contents) {
   YAML::Node OutStream        = YAML::Clone(Template);
   OutStream["UsePointerFile"] = false;
   OutStream["Filename"]       = OutFile;
   OutStream["IfExists"]       = "replace";
   OutStream["Precision"]      = "double";
   OutStream["Freq"]           = 1;
   OutStream["FreqUnits"]      = "OnDemand";
   OutStream["UseStartEnd"]    = false;
   OutStream.remove("PointerFilename");
   YAML::Node ContentsNode;
   ContentsNode.push_back(Contents);
   OutStream["Contents"] = ContentsNode;
   return OutStream;
}

//------------------------------------------------------------------------------
// Writes a copy of the input config in which HorzMeshIn reads MeshFile and
// adds a write stream for its contents. Only task 0 writes the file; all
// tasks then read the same copy. Returns the new config file name.
std::string makeRoundTripConfig(const std::string &InFile,
                                const std::string &MeshFile,
                                const std::string &HorzOutFile,
                                const std::string &ConfigFile, MachEnv *Env) {

   if (Env->getMyTask() == 0) {
      YAML::Node Root    = YAML::LoadFile(InFile);
      YAML::Node Streams = Root["Omega"]["IOStreams"];

      Streams["HorzMeshIn"]["Filename"] = MeshFile;
      Streams[HorzOutStreamName] =
          makeOutStream(Streams["RestartWrite"], HorzOutFile, "HorzMeshIn");

      std::ofstream Out(ConfigFile);
      Out << Root;
   }
   MPI_Barrier(Env->getComm());

   return ConfigFile;

} // end makeRoundTripConfig

//------------------------------------------------------------------------------
// Writes an output stream, which Omega validates (expanding group names into
// fields) before writing. Only the streams written here are validated because
// the fields of other streams (eg InitialState) are not defined in this test.
void writeStream(const std::string &StreamName, const Clock *ModelClock) {
   if (!IOStream::get(StreamName)->validate())
      ABORT_ERROR("MeshRoundTrip: unable to validate stream {}", StreamName);
   IOStream::write(StreamName, ModelClock, true);
}

//------------------------------------------------------------------------------
// The older MPAS name of a variable, which Omega also accepts when reading
// (see IOStream::readFieldData): the first letter is lower case and "Layer"
// is replaced by "Level"
std::string mpasName(const std::string &Name) {
   std::string OldName = Name;
   OldName[0]          = std::tolower(OldName[0]);
   size_t Pos          = OldName.find("Layer");
   if (Pos != std::string::npos)
      OldName.replace(Pos, 5, "Level");
   return OldName;
}

//------------------------------------------------------------------------------
// Basic information about a variable in a file
struct VarInfo {
   int VarID = -1;
   nc_type Type;
   std::vector<PIO_Offset> Shape;
   PIO_Offset Size   = 1;
   bool LeadingUnlim = false; // first dimension is the unlimited (time) dim
   std::string LeadingDimName;
};

VarInfo getVarInfo(int FileID, int VarID) {
   VarInfo Info;
   Info.VarID = VarID;
   int NDims;
   PIOc_inq_vartype(FileID, VarID, &Info.Type);
   PIOc_inq_varndims(FileID, VarID, &NDims);
   std::vector<int> DimIDs(NDims);
   if (NDims > 0)
      PIOc_inq_vardimid(FileID, VarID, DimIDs.data());
   for (int IDim = 0; IDim < NDims; ++IDim) {
      PIO_Offset Len;
      PIOc_inq_dimlen(FileID, DimIDs[IDim], &Len);
      Info.Shape.push_back(Len);
      Info.Size *= Len;
   }
   if (NDims > 0) {
      int UnlimDimID = -1;
      PIOc_inq_unlimdim(FileID, &UnlimDimID);
      Info.LeadingUnlim = DimIDs[0] == UnlimDimID;
      char NameBuf[PIO_MAX_NAME + 1];
      PIOc_inq_dimname(FileID, DimIDs[0], NameBuf);
      Info.LeadingDimName = NameBuf;
   }
   return Info;
}

std::string shapeString(const std::vector<PIO_Offset> &Shape) {
   std::string Str = "(";
   for (size_t I = 0; I < Shape.size(); ++I) {
      if (I > 0)
         Str += ", ";
      Str += std::to_string(Shape[I]);
   }
   return Str + ")";
}

bool isIntType(nc_type Type) {
   return Type == PIO_BYTE or Type == PIO_SHORT or Type == PIO_INT or
          Type == PIO_INT64;
}

bool isRealType(nc_type Type) {
   return Type == PIO_FLOAT or Type == PIO_DOUBLE;
}

//------------------------------------------------------------------------------
// Compares one output variable with the matching input variable and returns
// true if they are identical in shape and values. Integers are compared
// exactly as 64-bit integers and reals as doubles, except that an input
// variable written as single precision (for example, in a single-precision
// build) is first rounded to single precision.
//
// For connectivity, RowCounts holds the number of valid entries in each row
// (eg nEdgesOnCell for CellsOnCell) from the input file. Entries beyond it
// are padding, which mesh files fill with 0 or an invalid ID by different
// conventions, so they are not compared, only counted for information.
bool compareVar(const std::string &OutName, const std::string &InName,
                int OutFileID, const VarInfo &Out, int InFileID,
                const VarInfo &In,
                const std::vector<long long> *RowCounts = nullptr) {

   bool Identical = true;

   if (Out.Shape != In.Shape) {
      LOG_ERROR("MeshRoundTrip: {} has shape {} in the output but {} has "
                "shape {} in the input",
                OutName, shapeString(Out.Shape), InName, shapeString(In.Shape));
      Identical = false;

      // A static field written as time dependent has an extra leading
      // unlimited dimension of length 1 but the same values, so the values
      // can still be compared
      std::vector<PIO_Offset> Trailing;
      if (Out.LeadingUnlim and Out.Shape.size() == In.Shape.size() + 1 and
          Out.Shape[0] == 1)
         Trailing.assign(Out.Shape.begin() + 1, Out.Shape.end());
      if (!Trailing.empty() and Trailing == In.Shape) {
         LOG_ERROR("MeshRoundTrip: {} was written with the unlimited "
                   "dimension {} but is not time dependent in the input",
                   OutName, Out.LeadingDimName);
      } else {
         return Identical;
      }
   }

   bool IsInt  = isIntType(Out.Type) and isIntType(In.Type);
   bool IsReal = (isRealType(Out.Type) or isIntType(Out.Type)) and
                 (isRealType(In.Type) or isIntType(In.Type)) and !IsInt;
   if (!IsInt and !IsReal) {
      LOG_ERROR("MeshRoundTrip: {} (type {}) and {} (type {}) are not both "
                "numeric",
                OutName, Out.Type, InName, In.Type);
      return false;
   }
   if (Out.Type != In.Type) {
      LOG_INFO("MeshRoundTrip: {} has type {} in the output but {} has type "
               "{} in the input",
               OutName, Out.Type, InName, In.Type);
   }

   long Mismatches = 0;
   bool SameOffset = true; // all mismatches have the same output - input
   double Offset   = 0.0;

   // Whether entry I (of the flattened array) is padding beyond the count
   // of valid entries in its row
   PIO_Offset NCols = In.Shape.empty() ? 1 : In.Shape.back();
   if (RowCounts != nullptr and
       (In.Shape.size() != 2 or
        static_cast<PIO_Offset>(RowCounts->size()) != In.Shape[0])) {
      LOG_ERROR("MeshRoundTrip: row counts do not match the shape {} of {}",
                shapeString(In.Shape), InName);
      return false;
   }
   long PaddingDiffs = 0;
   auto isPadding    = [&](PIO_Offset I) {
      return RowCounts != nullptr and (I % NCols) >= (*RowCounts)[I / NCols];
   };

   auto record = [&](PIO_Offset I, double InV, double OutV) {
      if (Mismatches < MaxReport) {
         LOG_ERROR("MeshRoundTrip: {}[{}] input {} output {}", OutName, I, InV,
                   OutV);
      }
      if (Mismatches == 0)
         Offset = OutV - InV;
      else if (OutV - InV != Offset)
         SameOffset = false;
      ++Mismatches;
   };

   if (IsInt) {
      std::vector<long long> OutI(Out.Size);
      std::vector<long long> InI(In.Size);
      PIOc_get_var_longlong(OutFileID, Out.VarID, OutI.data());
      PIOc_get_var_longlong(InFileID, In.VarID, InI.data());
      for (PIO_Offset I = 0; I < Out.Size; ++I) {
         if (OutI[I] == InI[I])
            continue;
         if (isPadding(I))
            ++PaddingDiffs;
         else
            record(I, static_cast<double>(InI[I]),
                   static_cast<double>(OutI[I]));
      }
   } else {
      std::vector<double> OutVal(Out.Size);
      std::vector<double> InVal(In.Size);
      PIOc_get_var_double(OutFileID, Out.VarID, OutVal.data());
      PIOc_get_var_double(InFileID, In.VarID, InVal.data());
      bool Single = Out.Type == PIO_FLOAT;
      for (PIO_Offset I = 0; I < Out.Size; ++I) {
         double InV  = InVal[I];
         double OutV = OutVal[I];
         if (Single)
            InV = static_cast<double>(static_cast<float>(InV));
         bool BothNaN = std::isnan(InV) and std::isnan(OutV);
         if (InV != OutV and !BothNaN)
            record(I, InV, OutV);
      }
   }

   if (PaddingDiffs > 0) {
      LOG_INFO("MeshRoundTrip: {} differs from input {} in {} padding entries "
               "(beyond the count of valid entries), which are not compared",
               OutName, InName, PaddingDiffs);
   }

   if (Mismatches > 0) {
      LOG_ERROR("MeshRoundTrip: {} differs from input {} in {} of {} entries",
                OutName, InName, Mismatches, Out.Size);
      if (SameOffset) {
         LOG_ERROR("MeshRoundTrip: every mismatched entry of {} has output - "
                   "input = {}",
                   OutName, Offset);
      }
      Identical = false;
   } else if (Identical) {
      LOG_INFO("MeshRoundTrip: {} matches input {} ({} entries)", OutName,
               InName, Out.Size);
   } else {
      LOG_ERROR("MeshRoundTrip: values of {} match input {} ({} entries)",
                OutName, InName, Out.Size);
   }
   return Identical;

} // end compareVar

//------------------------------------------------------------------------------
// Reads the first of the given names (Omega, then MPAS) that is a global
// attribute of the file. Strings are converted to lower case, as HorzMesh
// does when reading them. Returns false if none of the names is found.
bool readStrAtt(int FileID, const std::vector<std::string> &Names,
                std::string &Val) {
   for (const auto &Name : Names) {
      if (IO::readMeta(Name, Val, FileID, IO::GlobalID).isSuccess()) {
         std::transform(Val.begin(), Val.end(), Val.begin(),
                        [](unsigned char C) { return std::tolower(C); });
         return true;
      }
   }
   return false;
}

bool readRealAtt(int FileID, const std::vector<std::string> &Names, R8 &Val) {
   for (const auto &Name : Names)
      if (IO::readMeta(Name, Val, FileID, IO::GlobalID).isSuccess())
         return true;
   return false;
}

//------------------------------------------------------------------------------
// Checks the global mesh attributes that HorzMesh reads (under either the
// Omega or the MPAS name). Every valid horizontal mesh file needs
// on_a_sphere and is_periodic, a mesh on a sphere also needs sphere_radius
// and a periodic planar mesh also needs x_period and y_period. A required
// attribute that is missing from the output, or any attribute whose value
// differs from the input, is an error. Returns the number of errors.
int checkMeshAttributes(int InFileID, int OutFileID) {

   const std::vector<std::string> SphereNames   = {"OnSphere", "on_a_sphere"};
   const std::vector<std::string> PeriodicNames = {"IsPeriodic", "is_periodic"};

   std::string InSphere;
   std::string InPeriodic;
   bool OnSphere =
       readStrAtt(InFileID, SphereNames, InSphere) and InSphere == "yes";
   bool Periodic = !OnSphere and
                   readStrAtt(InFileID, PeriodicNames, InPeriodic) and
                   InPeriodic == "yes";

   int NFail = 0;

   auto checkStr = [&](const std::vector<std::string> &Names, bool Required) {
      std::string InVal;
      std::string OutVal;
      bool InFound  = readStrAtt(InFileID, Names, InVal);
      bool OutFound = readStrAtt(OutFileID, Names, OutVal);
      if (!OutFound) {
         if (Required) {
            LOG_ERROR("MeshRoundTrip: required global attribute {} (or {}) "
                      "is missing from the output",
                      Names[1], Names[0]);
            ++NFail;
         }
      } else if (InFound and OutVal != InVal) {
         LOG_ERROR("MeshRoundTrip: global attribute {} is {} in the input "
                   "but {} in the output",
                   Names[1], InVal, OutVal);
         ++NFail;
      } else {
         LOG_INFO("MeshRoundTrip: global attribute {} = {} matches the input",
                  Names[1], OutVal);
      }
   };

   auto checkReal = [&](const std::vector<std::string> &Names, bool Required) {
      R8 InVal;
      R8 OutVal;
      bool InFound  = readRealAtt(InFileID, Names, InVal);
      bool OutFound = readRealAtt(OutFileID, Names, OutVal);
      if (!OutFound) {
         if (Required) {
            LOG_ERROR("MeshRoundTrip: required global attribute {} (or {}) "
                      "is missing from the output",
                      Names[1], Names[0]);
            ++NFail;
         }
      } else if (InFound and OutVal != InVal) {
         LOG_ERROR("MeshRoundTrip: global attribute {} is {} in the input "
                   "but {} in the output",
                   Names[1], InVal, OutVal);
         ++NFail;
      } else {
         LOG_INFO("MeshRoundTrip: global attribute {} = {} matches the input",
                  Names[1], OutVal);
      }
   };

   checkStr(SphereNames, true);
   checkStr(PeriodicNames, true);
   checkReal({"SphereRadius", "sphere_radius"}, OnSphere);
   checkReal({"XPeriod", "x_period"}, Periodic);
   checkReal({"YPeriod", "y_period"}, Periodic);

   return NFail;

} // end checkMeshAttributes

//------------------------------------------------------------------------------
// Reads an integer variable from the input file by its Omega name or else its
// MPAS name. Returns false if it is in neither.
bool readInputInts(int InFileID, const std::string &Name,
                   std::vector<long long> &Values) {
   int VarID;
   if (PIOc_inq_varid(InFileID, Name.c_str(), &VarID) != PIO_NOERR and
       PIOc_inq_varid(InFileID, mpasName(Name).c_str(), &VarID) != PIO_NOERR)
      return false;
   VarInfo Info = getVarInfo(InFileID, VarID);
   Values.resize(Info.Size);
   PIOc_get_var_longlong(InFileID, VarID, Values.data());
   return true;
}

//------------------------------------------------------------------------------
// Compares every variable in the output file with the matching variable in
// the input file. If CheckAttributes is true, also checks the global mesh
// attributes. Returns the number of variables and attributes that do not
// match.
int compareFiles(const std::string &InFile, const std::string &OutFile,
                 bool CheckAttributes) {

   LOG_INFO("MeshRoundTrip: comparing {} with {}", OutFile, InFile);

   int InFileID;
   int OutFileID;
   IO::openFileRead(InFileID, InFile);
   IO::openFileRead(OutFileID, OutFile);

   // Variables every Omega output stream writes that are not mesh variables
   const std::set<std::string> NonMeshVars = {"time", "time_bnds"};

   int NOutVars;
   PIOc_inq_nvars(OutFileID, &NOutVars);

   int NFail     = 0;
   int NCompared = 0;
   std::set<std::string> MatchedInVars;

   // The connectivity arrays whose rows are padded beyond a count of valid
   // entries, and the counts (from the input file) for each
   const std::map<std::string, std::string> RowCountNames = {
       {"CellsOnCell", "NEdgesOnCell"},
       {"EdgesOnCell", "NEdgesOnCell"},
       {"VerticesOnCell", "NEdgesOnCell"},
       {"EdgesOnEdge", "NEdgesOnEdge"},
       {"ReconStencilCell", "NEdgesReconOnCell"}};
   std::map<std::string, std::vector<long long>> RowCounts;
   for (const auto &Entry : RowCountNames) {
      const std::string &CountName = Entry.second;
      std::vector<long long> Counts;
      if (RowCounts.count(CountName) == 0 and
          readInputInts(InFileID, CountName, Counts))
         RowCounts[CountName] = Counts;
   }

   for (int OutVarID = 0; OutVarID < NOutVars; ++OutVarID) {
      char NameBuf[PIO_MAX_NAME + 1];
      PIOc_inq_varname(OutFileID, OutVarID, NameBuf, sizeof(NameBuf));
      std::string OutName = NameBuf;
      if (NonMeshVars.count(OutName) > 0)
         continue;

      // Find the input variable by the Omega name, then the MPAS name, just
      // as Omega does when reading
      std::string InName = OutName;
      int InVarID;
      if (PIOc_inq_varid(InFileID, InName.c_str(), &InVarID) != PIO_NOERR) {
         InName = mpasName(OutName);
         if (PIOc_inq_varid(InFileID, InName.c_str(), &InVarID) != PIO_NOERR) {
            // Not an error: Omega supplies defaults for some fields that
            // are missing from the mesh file and then writes them out
            LOG_WARN("MeshRoundTrip: output variable {} (or {}) is not in "
                     "the input file, so it is not compared",
                     OutName, InName);
            continue;
         }
      }
      MatchedInVars.insert(InName);

      VarInfo Out                          = getVarInfo(OutFileID, OutVarID);
      VarInfo In                           = getVarInfo(InFileID, InVarID);
      const std::vector<long long> *Counts = nullptr;
      auto CountName                       = RowCountNames.find(OutName);
      if (CountName != RowCountNames.end() and
          RowCounts.count(CountName->second) > 0)
         Counts = &RowCounts[CountName->second];
      bool Identical =
          compareVar(OutName, InName, OutFileID, Out, InFileID, In, Counts);
      ++NCompared;
      if (!Identical)
         ++NFail;
   }

   if (NCompared == 0) {
      LOG_ERROR("MeshRoundTrip: no mesh variables found in {}", OutFile);
      ++NFail;
   }

   // List the input variables Omega does not write (for information only,
   // since much of the mesh, such as the connectivity, is only used to build
   // the decomposition)
   int NInVars;
   PIOc_inq_nvars(InFileID, &NInVars);
   std::string NotWritten;
   for (int InVarID = 0; InVarID < NInVars; ++InVarID) {
      char NameBuf[PIO_MAX_NAME + 1];
      PIOc_inq_varname(InFileID, InVarID, NameBuf, sizeof(NameBuf));
      if (MatchedInVars.count(NameBuf) == 0)
         NotWritten += std::string(" ") + NameBuf;
   }
   LOG_INFO("MeshRoundTrip: compared {} variables; input variables not "
            "written:{}",
            NCompared, NotWritten);

   if (CheckAttributes)
      NFail += checkMeshAttributes(InFileID, OutFileID);

   IO::closeFile(OutFileID);
   IO::closeFile(InFileID);

   return NFail;

} // end compareFiles

} // namespace

//------------------------------------------------------------------------------
// The test driver for the mesh round trip
int main(int argc, char *argv[]) {

   std::string MeshFile = "OmegaMesh.nc";
   if (argc > 1)
      MeshFile = argv[1];

   // Distinguish the files of this test by mesh so that tests on different
   // meshes can run at the same time in the same directory
   std::string Stem        = std::filesystem::path(MeshFile).stem().string();
   std::string HorzOutFile = "MeshRoundTrip_" + Stem + "_horz.nc";
   std::string ConfigFile  = "omega.meshroundtrip." + Stem + ".yml";

   int NFail = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      MachEnv::init(MPI_COMM_WORLD);
      MachEnv *DefEnv  = MachEnv::getDefault();
      MPI_Comm DefComm = DefEnv->getComm();

      initLogging(DefEnv);
      LOG_INFO("------ Mesh round trip test for {} ------", MeshFile);

      Config("Omega");
      Config::readAll(makeRoundTripConfig("omega.yml", MeshFile, HorzOutFile,
                                          ConfigFile, DefEnv));

      // Read the mesh in the same order as Omega does in OceanInit
      TimeStepper::init1();
      Clock *ModelClock = TimeStepper::getDefault()->getClock();
      IO::init(DefComm);
      IOStream::init(ModelClock);
      Field::init(ModelClock);
      Decomp::init();
      Halo::init();
      HorzMesh::init(ModelClock);

      // Write the mesh back out through an Omega stream
      writeStream(HorzOutStreamName, ModelClock);
      MPI_Barrier(DefComm);

      // Compare the output with the mesh file
      NFail = compareFiles(MeshFile, HorzOutFile, true);

      IOStream::finalize();
      TimeStepper::clear();
      HorzMesh::clear();
      Field::clear();
      Dimension::clear();
      Halo::clear();
      Decomp::clear();
   }
   Pacer::finalize();
   Kokkos::finalize();

   if (NFail > 0)
      ABORT_ERROR(
          "------ Mesh round trip test FAIL: {} variables or attributes "
          "differ ------",
          NFail);
   LOG_INFO("------ Mesh round trip test successful ------");

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0;

} // end main
//===-----------------------------------------------------------------------===/
