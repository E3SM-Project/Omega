//===-- ocn/Eos.cpp - Equation of State ------------------*- C++ -*-===//
//
// The Eos class is responsible for managing the equation of state. It
// has a linear EOS and TEOS-10 EOS option, which is determined at
// initialization. It contains arrays that store the specific volume and
// displaced specific volume data.
//
//===----------------------------------------------------------------------===//

#include "Eos.h"
#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

/// Constructor for Teos10Eos
Teos10Eos::Teos10Eos(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

/// Constructor for LinearEos
LinearEos::LinearEos(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

/// Constructor for ConstantEos
ConstantEos::ConstantEos(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

/// Constructor for Teos10 squared Brunt-Vaisala frequency
Teos10BruntVaisalaFreqSq::Teos10BruntVaisalaFreqSq(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

/// Constructor for Linear squared Brunt-Vaisala frequency
LinearBruntVaisalaFreqSq::LinearBruntVaisalaFreqSq(const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      GeomZMid(VCoord->GeomZMid) {}

/// Constructor for Eos
Eos::Eos(const std::string &Name, ///< [in] Name for eos object
         const HorzMesh *Mesh,    ///< [in] Horizontal mesh
         const VertCoord *VCoord  ///< [in] Vertical coordinate
         )
    : ComputeSpecVolLinear(VCoord), ComputeSpecVolTeos10(VCoord),
      ComputeSpecVolConstant(VCoord), ComputeBruntVaisalaFreqSqLinear(VCoord),
      ComputeBruntVaisalaFreqSqTeos10(VCoord), Name(Name), Mesh(Mesh),
      VCoord(VCoord) {
   SpecVol = Array2DReal("SpecVol", Mesh->NCellsSize, VCoord->NVertLayers);
   SpecVolDisplaced =
       Array2DReal("SpecVolDisplaced", Mesh->NCellsSize, VCoord->NVertLayers);
   BruntVaisalaFreqSq = Array2DReal("BruntVaisalaFreqSq", Mesh->NCellsSize,
                                    VCoord->NVertLayersP1);
   SpecVolDCt =
       Array2DReal("SpecVolDCt", Mesh->NCellsSize, VCoord->NVertLayers);
   SpecVolDSa =
       Array2DReal("SpecVolDSa", Mesh->NCellsSize, VCoord->NVertLayers);
   SpecVolDP = Array2DReal("SpecVolDP", Mesh->NCellsSize, VCoord->NVertLayers);

   defineFields();
}

/// Destructor for Eos
Eos::~Eos() {}

/// Instance management
Eos *Eos::Instance = nullptr;

/// Get instance of Eos
Eos *Eos::getInstance() { return Instance; }

/// Destroy instance of Eos
void Eos::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initializes the Eos (Equation of State) class and its options.
/// it ASSUMES that HorzMesh was initialized and initializes the Eos class by
/// using the default mesh, reading the config file, and setting parameters
/// for either a Linear or TEOS-10 equation.
void Eos::init() {

   if (!Instance) {
      Instance =
          new Eos("Default", HorzMesh::getDefault(), VertCoord::getDefault());
   }

   Error Err; // error code

   /// Retrieve default eos
   Eos *eos = Eos::getInstance();

   /// Get EosConfig group from Omega config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config EosConfig("Eos");
   Err += OmegaConfig->get(EosConfig);
   CHECK_ERROR_ABORT(Err, "Eos::init: Eos group not found in Config");

   /// Get EosType from EosConfig
   /// and set the EosChoice accordingly
   std::string EosTypeStr;
   Err += EosConfig.get("EosType", EosTypeStr);
   CHECK_ERROR_ABORT(Err, "Eos::init: EosType subgroup not found in EosConfig");

   /// Set EosChoice and parameters based on EosTypeStr
   if (EosTypeStr == "Linear" or EosTypeStr == "linear") {
      Config EosLinConfig("Linear");
      Err += EosConfig.get(EosLinConfig);

      eos->EosChoice = EosType::LinearEos;

      CHECK_ERROR_ABORT(Err,
                        "Eos::init: Linear subgroup not found in EosConfig");
      Err += EosLinConfig.get("DRhoDT", eos->ComputeSpecVolLinear.DRhodT);
      CHECK_ERROR_ABORT(
          Err, "Eos::init: Parameter Linear:DRhodT not found in EosLinConfig");

      Err += EosLinConfig.get("DRhoDS", eos->ComputeSpecVolLinear.DRhodS);
      CHECK_ERROR_ABORT(
          Err, "Eos::init: Parameter Linear:DRhodS not found in EosLinConfig");

      Err += EosLinConfig.get("RhoT0S0", eos->ComputeSpecVolLinear.RhoT0S0);
      CHECK_ERROR_ABORT(
          Err, "Eos::init: Parameter Linear:RhoT0S0 not found in EosLinConfig");
   } else if ((EosTypeStr == "teos10") or (EosTypeStr == "teos-10") or
              (EosTypeStr == "TEOS-10")) {
      eos->EosChoice = EosType::Teos10Eos;
   } else if ((EosTypeStr == "constant") or (EosTypeStr == "Constant")) {
      eos->EosChoice = EosType::ConstantEos;
   } else {
      ABORT_ERROR("Eos::init: Unknown EosType requested");
   }
} // end init

/// Compute specific volume for all cells/layers (no displacement)
void Eos::computeSpecVol(const Array2DReal &ConservTemp,
                         const Array2DReal &AbsSalinity,
                         const Array2DReal &Pressure) {
   // count the evaluations this call performs, for the pressure gradient
   // cost check
   SpecVolEvalCount += static_cast<I8>(Mesh->NCellsAll) * VCoord->NVertLayers;

   OMEGA_SCOPE(LocSpecVol, SpecVol); /// Create a local view for computation
   OMEGA_SCOPE(LocComputeSpecVolLinear,
               ComputeSpecVolLinear); /// Local view for linear EOS computation
   OMEGA_SCOPE(LocComputeSpecVolTeos10,
               ComputeSpecVolTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(LocComputeSpecVolConstant,
               ComputeSpecVolConstant); /// Local view for constant computation
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   I4 KDisp = 0; /// No displacement in this case

   /// Dispatch to the correct EOS calculation
   if (EosChoice == EosType::LinearEos) {
      parallelForOuter(
          "eos-linear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolLinear(LocSpecVol, ICell, KChunk,
                                            ConservTemp, AbsSalinity);
                 });
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      parallelForOuter(
          "eos-teos10", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolTeos10(LocSpecVol, ICell, KChunk,
                                            ConservTemp, AbsSalinity, Pressure,
                                            KDisp);
                 });
          });
   } else if (EosChoice == EosType::ConstantEos) {
      parallelForOuter(
          "eos-constant", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolConstant(LocSpecVol, ICell, KChunk,
                                              ConservTemp, AbsSalinity);
                 });
          });
   }
}

/// Compute displaced specific volume (for vertical displacement)
void Eos::computeSpecVolDisp(const Array2DReal &ConservTemp,
                             const Array2DReal &AbsSalinity,
                             const Array2DReal &Pressure, I4 KDisp) {
   // count the evaluations this call performs, for the pressure gradient
   // cost check
   SpecVolEvalCount += static_cast<I8>(Mesh->NCellsAll) * VCoord->NVertLayers;

   OMEGA_SCOPE(LocSpecVolDisplaced,
               SpecVolDisplaced); /// Local view for computation
   OMEGA_SCOPE(LocComputeSpecVolLinear,
               ComputeSpecVolLinear); /// Local view for linear EOS computation
   OMEGA_SCOPE(LocComputeSpecVolTeos10,
               ComputeSpecVolTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(LocComputeSpecVolConstant,
               ComputeSpecVolConstant); /// Local view for constant computation
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   /// Dispatch to the correct EOS calculation
   /// If EosChoice is Linear, the displaced specific
   /// volume is the same as the specific volume
   if (EosChoice == EosType::LinearEos) {
      parallelForOuter(
          "eos-linear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolLinear(LocSpecVolDisplaced, ICell, KChunk,
                                            ConservTemp, AbsSalinity);
                 });
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      parallelForOuter(
          "eos-teos10", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolTeos10(LocSpecVolDisplaced, ICell, KChunk,
                                            ConservTemp, AbsSalinity, Pressure,
                                            KDisp);
                 });
          });
   } else if (EosChoice == EosType::ConstantEos) {
      parallelForOuter(
          "eos-constant", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolConstant(LocSpecVolDisplaced, ICell,
                                              KChunk, ConservTemp, AbsSalinity);
                 });
          });
   }
}

/// Compute specific volume and its first derivatives for all cells/layers
void Eos::computeSpecVolAndDerivs(const Array2DReal &ConservTemp,
                                  const Array2DReal &AbsSalinity,
                                  const Array2DReal &Pressure) {
   // count the evaluations this call performs, for the pressure gradient
   // cost check
   SpecVolEvalCount += static_cast<I8>(Mesh->NCellsAll) * VCoord->NVertLayers;

   OMEGA_SCOPE(LocSpecVol, SpecVol);       /// Local views for computation
   OMEGA_SCOPE(LocSpecVolDCt, SpecVolDCt); /// Temperature derivative
   OMEGA_SCOPE(LocSpecVolDSa, SpecVolDSa); /// Salinity derivative
   OMEGA_SCOPE(LocSpecVolDP, SpecVolDP);   /// Pressure derivative
   OMEGA_SCOPE(LocComputeSpecVolLinear,
               ComputeSpecVolLinear); /// Local view for linear EOS computation
   OMEGA_SCOPE(LocComputeSpecVolTeos10,
               ComputeSpecVolTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(LocComputeSpecVolConstant,
               ComputeSpecVolConstant); /// Local view for constant computation
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   /// Dispatch to the correct EOS calculation
   if (EosChoice == EosType::LinearEos) {
      parallelForOuter(
          "eos-derivs-linear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolLinear.calcSpecVolAndDerivsInChunk(
                        LocSpecVol, LocSpecVolDCt, LocSpecVolDSa, LocSpecVolDP,
                        ICell, KChunk, ConservTemp, AbsSalinity);
                 });
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      parallelForOuter(
          "eos-derivs-teos10", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolTeos10.calcSpecVolAndDerivsInChunk(
                        LocSpecVol, LocSpecVolDCt, LocSpecVolDSa, LocSpecVolDP,
                        ICell, KChunk, ConservTemp, AbsSalinity, Pressure);
                 });
          });
   } else if (EosChoice == EosType::ConstantEos) {
      parallelForOuter(
          "eos-derivs-constant", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeSpecVolConstant.calcSpecVolAndDerivsInChunk(
                        LocSpecVol, LocSpecVolDCt, LocSpecVolDSa, LocSpecVolDP,
                        ICell, KChunk, ConservTemp, AbsSalinity);
                 });
          });
   }
}

/// Compute squared Brunt-Vaisala frequency for all cells/layers
void Eos::computeBruntVaisalaFreqSq(const Array2DReal &ConservTemp,
                                    const Array2DReal &AbsSalinity,
                                    const Array2DReal &Pressure,
                                    const Array2DReal &SpecVol) {
   OMEGA_SCOPE(LocBruntVaisalaFreqSq,
               BruntVaisalaFreqSq); /// Local view for computation
   OMEGA_SCOPE(
       LocComputeBruntVaisalaFreqSqLinear,
       ComputeBruntVaisalaFreqSqLinear); /// Local view for linear computation
   OMEGA_SCOPE(
       LocComputeBruntVaisalaFreqSqTeos10,
       ComputeBruntVaisalaFreqSqTeos10); /// Local view for TEOS-10 computation
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   /// Dispatch to the correct squared Brunt-Vaisala frequency calculation
   if (EosChoice == EosType::LinearEos) {
      /// If Linear EOS, use linear squared Brunt-Vaisala frequency calculation
      parallelForOuter(
          "bvf-linear", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             // Compute Brunt-Vaisala frequency at interior vertical interfaces
             const int KMin   = MinLayerCell(ICell) + 1;
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeBruntVaisalaFreqSqLinear(LocBruntVaisalaFreqSq,
                                                       ICell, KChunk, SpecVol);
                 });

             teamBarrier(Team);

             // Fill Brunt-Vaisala frequency at vertical boundaries using the
             // closest valid value. This is equivalent to doing one-sided
             // differencing at the boundary.
             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    LocBruntVaisalaFreqSq(ICell, MinLayerCell(ICell)) =
                        LocBruntVaisalaFreqSq(ICell, KMin);
                    LocBruntVaisalaFreqSq(ICell, MaxLayerCell(ICell) + 1) =
                        LocBruntVaisalaFreqSq(ICell, KMax);
                 });
          });
   } else if (EosChoice == EosType::Teos10Eos) {
      /// If TEOS-10 EOS, use TEOS-10 squared Brunt-Vaisala frequency
      /// calculation
      parallelForOuter(
          "bvf-teos10", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
             // Compute Brunt-Vaisala frequency at interior vertical interfaces
             const int KMin   = MinLayerCell(ICell) + 1;
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocComputeBruntVaisalaFreqSqTeos10(
                        LocBruntVaisalaFreqSq, ICell, KChunk, ConservTemp,
                        AbsSalinity, Pressure, SpecVol);
                 });

             teamBarrier(Team);

             // Fill Brunt-Vaisala frequency at vertical boundaries using the
             // closest valid value. This is equivalent to doing one-sided
             // differencing at the boundary.
             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    LocBruntVaisalaFreqSq(ICell, MinLayerCell(ICell)) =
                        LocBruntVaisalaFreqSq(ICell, KMin);
                    LocBruntVaisalaFreqSq(ICell, MaxLayerCell(ICell) + 1) =
                        LocBruntVaisalaFreqSq(ICell, KMax);
                 });
          });
   }
}

Real Eos::calcPtFromCt(const Real &Sa, const Real &Ct) const {
   if (EosChoice == EosType::Teos10Eos) {
      return ComputeSpecVolTeos10.calcPtFromCt(Sa, Ct);
   }

   return Ct;
}

Real Eos::calcCtFromPt(const Real &Sa, const Real &Pt) const {
   if (EosChoice == EosType::Teos10Eos) {
      return ComputeSpecVolTeos10.calcCtFromPt(Sa, Pt);
   }

   return Pt;
}

/// Define IO fields and metadata for output
void Eos::defineFields() {

   /// Set field names (append Name if not default)
   SpecVolFldName            = "SpecVol";
   SpecVolDisplacedFldName   = "SpecVolDisplaced";
   BruntVaisalaFreqSqFldName = "BruntVaisalaFreqSq";
   SpecVolDCtFldName         = "SpecVolDCt";
   SpecVolDSaFldName         = "SpecVolDSa";
   SpecVolDPFldName          = "SpecVolDP";
   if (Name != "Default") {
      SpecVolFldName.append(Name);
      SpecVolDisplacedFldName.append(Name);
      BruntVaisalaFreqSqFldName.append(Name);
      SpecVolDCtFldName.append(Name);
      SpecVolDSaFldName.append(Name);
      SpecVolDPFldName.append(Name);
   }

   /// Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLayers";

   /// Create and register the specific volume field
   auto SpecVolField =
       Field::create(SpecVolFldName,                   // Field name
                     "Layer-averaged Specific Volume", // Long Name
                     "m3 kg-1",                        // Units
                     "sea_water_specific_volume",      // CF-ish Name
                     0.0,                              // Min valid value
                     std::numeric_limits<Real>::max(), // Max valid value
                     NDims,                            // Number of dimensions
                     DimNames                          // Dimension names
       );
   /// Create and register the displaced specific volume field
   auto SpecVolDisplacedField =
       Field::create(SpecVolDisplacedFldName, // Field name
                     "Specific Volume displaced adiabatically "
                     "to specified layer",                  // long Name
                     "m3 kg-1",                             // Units
                     "sea_water_specific_volume_displaced", // CF-ish Name
                     0.0,                                   // Min valid value
                     std::numeric_limits<Real>::max(),      // Max valid value
                     NDims,   // Number of dimensions
                     DimNames // Dimension names
       );

   /// The specific volume derivatives are legitimately negative, so their
   /// valid range spans the full range of Real rather than starting at zero
   auto SpecVolDCtField = Field::create(
       SpecVolDCtFldName, // Field name
       "Derivative of specific volume with respect to conservative "
       "temperature",    // Long Name
       "m3 kg-1 degC-1", // Units
       // CF-ish Name
       "sea_water_specific_volume_derivative_wrt_conservative_temperature",
       std::numeric_limits<Real>::lowest(), // Min valid value
       std::numeric_limits<Real>::max(),    // Max valid value
       NDims,                               // Number of dimensions
       DimNames                             // Dimension names
   );

   auto SpecVolDSaField = Field::create(
       SpecVolDSaFldName, // Field name
       "Derivative of specific volume with respect to absolute "
       "salinity", // Long Name
       "m3 g-1",   // Units
       // CF-ish Name
       "sea_water_specific_volume_derivative_wrt_absolute_salinity",
       std::numeric_limits<Real>::lowest(), // Min valid value
       std::numeric_limits<Real>::max(),    // Max valid value
       NDims,                               // Number of dimensions
       DimNames                             // Dimension names
   );

   auto SpecVolDPField = Field::create(
       SpecVolDPFldName,                                         // Field name
       "Derivative of specific volume with respect to pressure", // Long Name
       "m3 kg-1 Pa-1",                                           // Units
       "sea_water_specific_volume_derivative_wrt_pressure",      // CF-ish Name
       std::numeric_limits<Real>::lowest(), // Min valid value
       std::numeric_limits<Real>::max(),    // Max valid value
       NDims,                               // Num dimensions
       DimNames                             // Dimension names
   );

   // Brunt-Vaisala frequency is located at interfaces
   DimNames[1] = "NVertLayersP1";

   /// Create and register the BruntVaisalaFreqSq field
   auto BruntVaisalaFreqSqField =
       Field::create(BruntVaisalaFreqSqFldName,                   // Field name
                     "Brunt-Vaisala frequency squared",           // Long Name
                     "s-2",                                       // Units
                     "sea_water_brunt_vaisala_frequency_squared", // CF-ish Name
                     std::numeric_limits<Real>::min(), // Min valid value
                     std::numeric_limits<Real>::max(), // Max valid value
                     NDims,                            // Number of dimensions
                     DimNames                          // Dimension names
       );

   // Create a field group for the eos-specific state fields
   EosGroupName = "Eos";
   if (Name != "Default") {
      EosGroupName.append(Name);
   }
   auto EosGroup = FieldGroup::create(EosGroupName);

   // Add fields to the EOS group
   EosGroup->addField(SpecVolDisplacedFldName);
   EosGroup->addField(SpecVolFldName);
   EosGroup->addField(BruntVaisalaFreqSqFldName);
   EosGroup->addField(SpecVolDCtFldName);
   EosGroup->addField(SpecVolDSaFldName);
   EosGroup->addField(SpecVolDPFldName);

   // Attach Kokkos views to the fields
   SpecVolDisplacedField->attachData<Array2DReal>(SpecVolDisplaced);
   SpecVolField->attachData<Array2DReal>(SpecVol);
   BruntVaisalaFreqSqField->attachData<Array2DReal>(BruntVaisalaFreqSq);
   SpecVolDCtField->attachData<Array2DReal>(SpecVolDCt);
   SpecVolDSaField->attachData<Array2DReal>(SpecVolDSa);
   SpecVolDPField->attachData<Array2DReal>(SpecVolDP);

} // end defineIOFields

} // namespace OMEGA
