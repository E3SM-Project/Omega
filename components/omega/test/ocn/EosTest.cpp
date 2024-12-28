//===-- Test driver for OMEGA GSW-C library -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA GSW-C external library
///
/// This driver tests that the GSW-C library can be called
/// and returns expected value (as published in Roquet et al 2015)
//
//===-----------------------------------------------------------------------===/

#include "Eos.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "EosConstants.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "mpi.h"

// added for debug
#include "AuxiliaryState.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"

#include <gswteos-10.h>

using namespace OMEGA;

constexpr Geometry Geom   = Geometry::Spherical;
constexpr int NVertLevels = 60;
// published values (TEOS-10) to test against
const Real TeosExpValueDelta = 0.0009776149797;
const Real TeosExpValueVol   = 0.0009732819628;
double Sa                    = 30.;
double Ct                    = 10.;
double P                     = 1000.;
const Real RTol              = 1e-10;
// expected value for Linear eos (and default parameters)
const Real LinearExpValue = 0.0009784735812133072;

//------------------------------------------------------------------------------
// The initialization routine for Eos testing. It calls various
// init routines, including the creation of the default decomposition.
I4 initEosTest(const std::string &mesh) {

   I4 Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("… in initEosTest");

   // Open config file
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_ERROR("Eos: Error reading config file");
      return Err;
   }

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing default decomposition");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing default mesh");
   }

   const auto &Mesh = HorzMesh::getDefault();
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);
   return Err;
}

int testEosMapping() {
   int Err = 0;
   // test initialization
   int EosErr = Eos::init();
   if (EosErr != 0) {
      Err++;
      LOG_ERROR("EosTest: error initializing default Eos");
   }

   // test retrievel of default
   Eos *DefEos = Eos::getDefault();

   if (DefEos) {
      LOG_INFO("EosTest: Default Eos retrieval PASS");
   } else {
      Err++;
      LOG_INFO("EosTest: Default Eos retrieval FAIL");
      return -1;
   }

   const auto *Mesh = HorzMesh::getDefault();
   // test creation of another Eos
   Eos::create("TestEos", Mesh, NVertLevels);

   if (Eos::get("TestEos")) {
      LOG_INFO("EosTest: Non-default Eos retrieval PASS");
   } else {
      Err++;
      LOG_INFO("EosTest: Non-default Eos retrieval FAIL");
   }

   // test erase
   Eos::erase("TestEos");

   if (Eos::get("TestEos")) {
      Err++;
      LOG_INFO("EosTest: Non-default Eos erase FAIL");
   } else {
      LOG_INFO("EosTest: Non-default Eos erase PASS");
   }

   // Eos::clear();

   return Err;
}

int testEosLinear() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   // create Eos to test
   Eos::create("LinearEos", Mesh, NVertLevels);
   Eos *TestEos       = Eos::get("LinearEos");
   TestEos->EosChoice = EosType::Linear;

   // create ocean state array
   Array2DReal Sarray = Array2DReal("Sarray", Mesh->NCellsAll, NVertLevels);
   Array2DReal Tarray = Array2DReal("Tarray", Mesh->NCellsAll, NVertLevels);
   Array2DReal Parray = Array2DReal("Parray", Mesh->NCellsAll, NVertLevels);
   // Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(Sarray, Sa);
   deepCopy(Tarray, Ct);
   deepCopy(Parray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   // Key calculation
   TestEos->computeSpecVol(TestEos->SpecVol, Tarray, Sarray, Parray);
   // LOG_INFO("EosTest, testEosLinear: produced SpecVol from Eos");

   // check on all array values
   int numMismatches = 0;
   Kokkos::parallel_reduce(
       "CheckSpecVolMatrix-linear",
       Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                              {Mesh->NCellsAll, NVertLevels}),
       KOKKOS_LAMBDA(const int i, const int j, int &localCount) {
          if (!isApprox(TestEos->SpecVol(i, j), LinearExpValue, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   bool allMatch = (numMismatches == 0);
   if (!allMatch) {
      Err++;
      LOG_ERROR("EosTest: SpecVol Linear isApprox FAIL, "
                "expected {}, got {} mismatches",
                LinearExpValue, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalc Linear: PASS");
   }

   Eos::erase("LinearEos");
   return Err;
}

int testEosTeos10() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   // create Eos to test
   Eos::create("TeosEos", Mesh, NVertLevels);
   Eos *TestEos       = Eos::get("TeosEos");
   TestEos->EosChoice = EosType::TEOS10Poly75t;

   // create ocean state array
   Array2DReal Sarray = Array2DReal("Sarray", Mesh->NCellsAll, NVertLevels);
   Array2DReal Tarray = Array2DReal("Tarray", Mesh->NCellsAll, NVertLevels);
   Array2DReal Parray = Array2DReal("Parray", Mesh->NCellsAll, NVertLevels);
   // Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(Sarray, Sa);
   deepCopy(Tarray, Ct);
   deepCopy(Parray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   // Key calculation
   TestEos->computeSpecVol(TestEos->SpecVol, Tarray, Sarray, Parray);
   // LOG_INFO("EosTest - testEosTeos: produced SpecVol from Eos");

   // check on all array values
   int numMismatches = 0;
   // LOG_INFO("Loop extent: {}, {}", Mesh->NCellsAll, NVertLevels);
   Kokkos::parallel_reduce(
       "CheckSpecVolMatrix-teos",
       Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                              {Mesh->NCellsAll, NVertLevels}),
       KOKKOS_LAMBDA(const int i, const int j, int &localCount) {
          if (!isApprox(TestEos->SpecVol(i, j), TeosExpValueVol, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   bool allMatch = (numMismatches == 0);
   if (!allMatch) {
      Err++;
      LOG_ERROR("EosTest: Teos SpecVol isApprox FAIL, "
                "expected {}, got {} mismatches",
                TeosExpValueVol, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalc TEOS-10: PASS");
   }

   Eos::erase("TeosEos");
   return Err;
}

void finalizeEosTest() {
   Eos::clear();
   HorzMesh::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

int gswcSpecVolCheckValue() {
   int Err         = 0;
   const Real RTol = 1e-10;

   double SpecVol = gsw_specvol(Sa, Ct, P);
   //  LOG_INFO("gswcSpecVolCheckValue: produced SpecVol from GSW-C module");
   //  LOG_INFO("Value of SpecVol: {}", SpecVol);
   bool Check = isApprox(SpecVol, TeosExpValueVol, RTol);
   if (!Check) {
      Err++;
      LOG_ERROR(
          "gswcSpecVolCheckValue: SpecVol isApprox FAIL, expected {}, got {}",
          TeosExpValueVol, SpecVol);
   }
   if (Err == 0) {
      LOG_INFO("gswcSpecVolCheckValue: PASS");
   }
   return Err;
}

// intermediate test accessing the coefficients
int fetchCoeff() {
   int Err       = 0;
   double ExpVal = 0.0010769995862;
   GSW_SPECVOL_COEFFICIENTS;
   // LOG_INFO("EosTest: called GSW_SPECVOL_COEFFICIENTS");
   // LOG_INFO("Value of V000: {}", V000);
   if (!isApprox(V000, ExpVal, RTol)) {
      Err++;
      LOG_ERROR("EosTest: Coeff V000 isApprox FAIL, expected {}, got {}",
                ExpVal, V000);
   }
   if (Err == 0) {
      LOG_INFO("GswcTeosTest: check PASS");
   }
   return Err;
}

int poly75tDeltaCheckValue() {
   int Err     = 0;
   const int K = 0;

   TEOS10Poly75t specvolpoly75t(NVertLevels);
   specvolpoly75t.calcPCoeffs(specvolpoly75t.specVolPcoeffs, K, Ct, Sa);
   Real Delta = specvolpoly75t.calcDelta(specvolpoly75t.specVolPcoeffs, K, P);
   // LOG_INFO("Teos10 poly75tDeltaCheckValue: produced delta from poly75t");
   // LOG_INFO("Value of Delta: {}", Delta);
   bool Check = isApprox(Delta, TeosExpValueDelta, RTol);
   if (!Check) {
      Err++;
      LOG_ERROR("Teos10 poly75tDeltaCheckValue: Delta isApprox FAIL, expected "
                "{}, got {}",
                TeosExpValueDelta, Delta);
   }
   if (Err == 0) {
      LOG_INFO("Teos10 poly75tDeltaCheckValue: PASS");
   }
   return Err;
}

int poly75tSpecVolCheckValue() {
   int Err             = 0;
   const Real RTol     = 1e-10;
   Array2DReal Sarray  = Array2DReal("Sarray", 1, 1);
   Array2DReal Tarray  = Array2DReal("Tarray", 1, 1);
   Array2DReal Parray  = Array2DReal("Parray", 1, 1);
   const int K         = 0;
   const int ICell     = 0;
   Array2DReal SpecVol = Array2DReal("SpecVol", 1, 1);

   Sarray(0, 0)  = Sa;
   Tarray(0, 0)  = Ct;
   Parray(0, 0)  = P;
   SpecVol(0, 0) = 0.0;

   TEOS10Poly75t specvolpoly75t(NVertLevels);
   specvolpoly75t(SpecVol, ICell, K, Tarray, Sarray, Parray);
   // LOG_INFO("Teos10 poly75tSpecVolCheckValue: produced SpecVol from
   // poly75t"); LOG_INFO("Value of SpecVol: {}", SpecVol(0, 0));
   bool Check = isApprox(SpecVol(0, 0), TeosExpValueVol, RTol);
   if (!Check) {
      Err++;
      LOG_ERROR("Teos10 poly75tSpecVolCheckValue: SpecVol isApprox FAIL, "
                "expected {}, got {}",
                TeosExpValueVol, SpecVol(0, 0));
   }
   if (Err == 0) {
      LOG_INFO("Teos10 poly75tSpecVolCheckValue: PASS");
   }
   return Err;
}

int linearSpecVolCheckValue() {
   int Err             = 0;
   const Real RTol2    = 1e-2; // (>> machine precision)
   Array2DReal Sarray  = Array2DReal("Sarray", 1, 1);
   Array2DReal Tarray  = Array2DReal("Tarray", 1, 1);
   const int K         = 0;
   const int ICell     = 0;
   Array2DReal SpecVol = Array2DReal("SpecVol", 1, 1);

   Sarray(0, 0)  = Sa;
   Tarray(0, 0)  = Ct;
   SpecVol(0, 0) = 0.0;

   LinearEOS specvollinear;
   specvollinear(SpecVol, ICell, K, Tarray, Sarray);
   // LOG_INFO("linearSpecVolCheckValue: produced SpecVol from linear EOS");
   // LOG_INFO("Value of SpecVol: {}", SpecVol(0, 0));
   bool Check = isApprox(SpecVol(0, 0), TeosExpValueVol, RTol); // expect False
   bool CheckClose =
       isApprox(SpecVol(0, 0), TeosExpValueVol, RTol2); // expect True
   if (Check) {
      Err++;
      LOG_ERROR("linearSpecVolCheckValue: SpecVol Linear is undistinguishable "
                "from TEOS10 Ref Value");
   }
   if (!CheckClose) {
      Err++;
      LOG_ERROR("linearSpecVolCheckValue: SpecVol TEOS10 and Linear are NOT "
                "close. Check input values");
   }
   if (Err == 0) {
      LOG_INFO("linearSpecVolCheckValue: PASS");
   }
   return Err;
}

int linearDensityLinearityTest() {
   int Err              = 0;
   const Real RTol      = 1e-10;
   const Real ExpDT     = 1.0;
   const Real ExpDS     = 2.4;
   Array2DReal Sarray1  = Array2DReal("Sarray1", 1, 1);
   Array2DReal Tarray1  = Array2DReal("Tarray1", 1, 1);
   Array2DReal Sarray2  = Array2DReal("Sarray2", 1, 1);
   Array2DReal Tarray2  = Array2DReal("Tarray2", 1, 1);
   const int K          = 0;
   const int ICell      = 0;
   Array2DReal SpecVol1 = Array2DReal("SpecVol1", 1, 1);
   Array2DReal SpecVol2 = Array2DReal("SpecVol2", 1, 1);
   Array2DReal SpecVol3 = Array2DReal("SpecVol3", 1, 1);
   Array2DReal SpecVol4 = Array2DReal("SpecVol4", 1, 1);

   Sarray1(0, 0)  = Sa;
   Tarray1(0, 0)  = 15.;
   Sarray2(0, 0)  = 33.;
   Tarray2(0, 0)  = Ct;
   SpecVol1(0, 0) = 0.0;
   SpecVol2(0, 0) = 0.0;
   SpecVol3(0, 0) = 0.0;
   SpecVol4(0, 0) = 0.0;

   LinearEOS specvollinear;
   specvollinear(SpecVol1, ICell, K, Tarray1, Sarray1);
   specvollinear(SpecVol2, ICell, K, Tarray2, Sarray1);
   specvollinear(SpecVol3, ICell, K, Tarray1, Sarray2);
   specvollinear(SpecVol4, ICell, K, Tarray2, Sarray2);
   Real DRhoT  = 1. / SpecVol2(0, 0) - 1. / SpecVol1(0, 0);
   Real DRhoS  = 1. / SpecVol3(0, 0) - 1. / SpecVol1(0, 0);
   Real DRhoTS = 1. / SpecVol4(0, 0) - 1. / SpecVol1(0, 0);
   // LOG_INFO("linearDensityLinearityTest: produced SpecVol from linear EOS");
   bool Check1 = isApprox(DRhoS, ExpDS, RTol);
   if (!Check1) {
      Err++;
      LOG_ERROR("linearDensityLinearityTest: DRhoS {}; expected {}", DRhoS,
                ExpDS);
   }
   bool Check2 = isApprox(DRhoT, ExpDT, RTol);
   if (!Check2) {
      Err++;
      LOG_ERROR("linearDensityLinearityTest: DRhoT {}; expected {}", DRhoT,
                ExpDT);
   }
   bool Check3 = isApprox(DRhoTS, DRhoS + DRhoT, RTol);
   if (!Check3) {
      Err++;
      LOG_ERROR("linearDensityLinearityTest: Sum(DRho) {}, DRhoTS {}",
                DRhoS + DRhoT, DRhoTS);
   }
   // LOG_INFO("linearDensityLinearityTest: Sum(DRho) is undistinguishable from"
   //          "DRhoTS");
   if (Err == 0) {
      LOG_INFO("linearDensityLinearityTest: PASS");
   }
   return Err;
}

// the main test (all in one to have the same log)
// Single value tests:
// --> one test calls the external GSW-C library
// and compares the specific volume to the published value
// --> next tests call the local Poly75t TEOS-10 calc
// --> next tests call the linear Eos
// Full array tests:
// --> one tests the initialization/retrieval of Eos
// --> next checks the value on a Eos with linear option
// --> next checks the value on a Eos with Teos-10 option
int eosTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initEosTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("EosTest: Error initializing");
   }
   const auto &Mesh = HorzMesh::getDefault();
   LOG_INFO("Single value checks:");
   Err += gswcSpecVolCheckValue();
   Err += fetchCoeff();
   // Err += poly75tDeltaCheckValue();
   // Err += poly75tSpecVolCheckValue();
   // Err += linearSpecVolCheckValue();
   // Err += linearDensityLinearityTest();
   LOG_INFO("Full array checks:");
   Err += testEosMapping();
   Err += testEosLinear();
   Err += testEosTeos10();

   if (Err == 0) {
      LOG_INFO("EosTest: Successful completion");
   }
   finalizeEosTest();

   return Err;
}

//------------------------------------------------------------------------------
// The test driver for Eos testing
int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   { RetVal += eosTest(); }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
