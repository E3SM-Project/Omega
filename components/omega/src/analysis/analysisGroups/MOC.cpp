//===-- analysis/analysisGroups/MOC.cpp -------------------------*- C++ -*-===//
//
// Implementation of MOC constructor and helper methods. Reads configuration
// for binning parameters, regions, and output specifications. Builds operator
// chains for each region to compute the MOC streamfunction, stores metadata
// for stream creation, and invokes base class method to create IOStreams.
//
//===----------------------------------------------------------------------===//

#include "analysisGroups/MOC.h"
#include "AnalysisOperator.h"
#include <iostream>

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructs a MOC analysis group by reading configuration, building operator
// chains for each region and transect, and creating IOStreams.
// Latitude-binned MOC streamfunction is computed by:
//   1. Converting VerticalPseudoVelocity to geometric coords
//   (PseudoToGeometric)
//   2. Computing vertical flux (BinaryMultiply with AreaCell)
//   3. Optionally masking to a region (ExtractRegion)
//   4. Accumulating into latitude bins (BinnedAccumulator)
//   5. Integrating horizontally south-to-north (PrefixSum)
//   6. Converting to Sverdrups (ScalarMultiply ×1e-6)
// Transect-based MOC streamfunction is computed by:
//   1. Converting PseudoThickness to geometric layer thickness
//   (PseudoToGeometric)
//   2. Computing edge transport (BinaryMultiply with NormalVelocity, then
//   DvEdge)
//   3. Accumulating across transect edges (TransectAccumulator)
//   4. Integrating vertically bottom-to-top (PrefixSum, reverse)
//   5. Converting to Sverdrups (ScalarMultiply ×1e-6)
MOC::MOC(const std::string &GroupName, Config &AnalysisGroupOptions,
         Analysis *AnalysisManager) {

   Error Err;

   // Read binning parameters with defaults
   I4 NumBins = 180; // Default: 1-degree bins
   Err        = AnalysisGroupOptions.get("NumBins", NumBins);
   if (Err.isFail()) {
      LOG_INFO("MOC: NumBins not specified, using default: {}", NumBins);
   }

   Real MinLat = -90.0; // Default: full latitude range
   Err         = AnalysisGroupOptions.get("MinLat", MinLat);
   if (Err.isFail()) {
      LOG_INFO("MOC: MinLat not specified, using default: {}", MinLat);
   }

   Real MaxLat = 90.0;
   Err         = AnalysisGroupOptions.get("MaxLat", MaxLat);
   if (Err.isFail()) {
      LOG_INFO("MOC: MaxLat not specified, using default: {}", MaxLat);
   }

   // Convert latitude from degrees to radians
   const Real DegToRad = M_PI / 180.0;
   Real MinLatRad      = MinLat * DegToRad;
   Real MaxLatRad      = MaxLat * DegToRad;

   // Read optional region list for regional MOC
   std::vector<std::string> RegionList;
   Err = AnalysisGroupOptions.get("Regions", RegionList);
   if (Err.isFail()) {
      // Default: compute global MOC only
      RegionList = {"Global"};
      LOG_INFO("MOC: Regions not specified, computing Global MOC only");
   }

   // Read optional transect list for transect-based MOC
   std::vector<std::string> TransectList;
   Err = AnalysisGroupOptions.get("Transects", TransectList);
   if (Err.isFail()) {
      // Default: no transect MOC
      TransectList = {};
      LOG_INFO(
          "MOC: Transects not specified, no transect MOC will be computed");
   }

   parseTemporalPeriods(AnalysisGroupOptions);

   // Build operator chains for each region
   std::vector<std::string> ChainStems;
   for (const auto &RegionName : RegionList) {

      // Determine mask field name for regional MOC
      // For global MOC, use empty string (no mask)
      std::string RegionMaskName;
      if (RegionName != "Global") {
         RegionMaskName = "RegionMask" + RegionName;
      }

      // Build the MOC operator chain for this region
      // This creates the chain string and stores config in ChainConfigs vector
      std::string ChainStem =
          buildMOCChain(RegionName, RegionMaskName, NumBins, MinLatRad,
                        MaxLatRad, AnalysisManager);
      ChainStems.push_back(ChainStem);
   }

   // Build operator chains for each transect
   for (const auto &TransectName : TransectList) {

      // Build the transect MOC operator chain
      // Mask and sign field names are derived from TransectName by
      // Analysis::parseChainAndBuildOps when it processes the
      // TransectAccumulator operator token.
      // This creates the chain string and stores config in ChainConfigs vector
      std::string ChainStem =
          buildTransectMOCChain(TransectName, AnalysisManager);
      ChainStems.push_back(ChainStem);
   }

   // Build temporal chains with MOC-specific operator configurations
   buildTemporalChains(ChainStems, AnalysisGroupOptions, AnalysisManager,
                       ChainConfigs);

   // Create IOStreams organized by output frequency and associate operators
   createAnalysisGroupStreams(GroupName, AnalysisGroupOptions, AnalysisManager);

} // end MOC constructor

//------------------------------------------------------------------------------
// Builds the complete operator chain string for computing MOC streamfunction
// for a single region. The chain consists of:
//   Initialization (once):
//     0. LatCell_CoordinateBinning: Assigns cells to latitude bins (static)
//   Runtime (each timestep):
//     1. VerticalPseudoVelocity_PseudoToGeometric: Convert to geometric coords
//     2. VerticalVelocity_BinaryMultiply(AreaCell): Compute vertical flux
//     3. ExtractRegion(RegionName): Apply regional mask (regional MOC only)
//     4. VerticalFlux_BinnedAccumulator(LatCell_BinIndex): Accumulate into bins
//     5. BinnedFlux_PrefixSum: Horizontal integration (south to north)
//     6. MOC_ScalarMultiply: Convert to Sverdrups (Sv)
//
// This method constructs the complete chain string that will later be parsed
// by buildTemporalChains(), which appends temporal operators and calls
// parseChainAndBuildOps().
//
// Returns the complete chain string (stem) for temporal operator appending
std::string MOC::buildMOCChain(const std::string &RegionName,
                               const std::string &RegionMaskName, I4 NumBins,
                               Real MinLat, Real MaxLat,
                               Analysis *AnalysisManager) {

   // -------------------------------------------------------------------------
   // Build the complete operator chain string
   // -------------------------------------------------------------------------
   // The complete chain encodes all transformation steps:
   //
   // For Global MOC:
   // VerticalPseudoVelocity_PseudoToGeometric_BinaryMultiply(AreaCell)_
   // BinnedAccumulator(LatCell_BinIndex)_PrefixSum_ScalarMultiply(1e-6)
   //
   // For Regional MOC:
   // VerticalPseudoVelocity_PseudoToGeometric_BinaryMultiply(AreaCell)_
   // ExtractRegion(RegionMask)_BinnedAccumulator(LatCell_BinIndex)_
   // PrefixSum_ScalarMultiply(1e-6)
   //
   // This will be parsed by Analysis::parseChainAndBuildOps() which:
   // 1. Splits by underscores to extract field names and operator types
   // 2. Creates operator instances with appropriate configurations
   // 3. Chains them together in sequence

   // Start with base vertical velocity field and conversion to geometric coords
   std::string ChainStr = "VerticalPseudoVelocity_PseudoToGeometric";

   // Multiply by AreaCell to get flux
   ChainStr += "_BinaryMultiply(AreaCell)";

   // Insert ExtractRegionOp for regional MOC after flux computation
   if (!RegionMaskName.empty()) {
      ChainStr += "_ExtractRegion(" + RegionMaskName + ")";
   }

   // Accumulate into latitude bins (uses LatCell_BinIndex from initialization)
   ChainStr += "_BinnedAccumulator(LatCell_BinIndex)";

   // Horizontal integration (south to north across latitude bins)
   ChainStr += "_PrefixSum";

   // Convert to Sverdrups
   ChainStr += "_ScalarMultiply(1.0e-6)";

   // -------------------------------------------------------------------------
   // Store configuration for this chain
   // -------------------------------------------------------------------------
   // Each operator in the chain needs configuration. Build a merged config
   // containing parameters for all operators using makeOpConfig/opParam
   // helpers.

   // Build config with all operator parameters
   Config ChainConfig = makeOpConfig(
       opParam("NumBins", NumBins), // BinnedAccumulator: number of bins
       opParam(
           "Dimension",
           0), // PrefixSum: integrate along bin/latitude dimension (horizontal)
       opParam("Reverse",
               false), // PrefixSum: forward integration (south-to-north)
       opParam("IOName", std::string("MOC_streamfunction_") + RegionName));

   // Append config to vector (indexed in same order as RegionList/ChainStems)
   ChainConfigs.push_back(ChainConfig);

   // -------------------------------------------------------------------------
   // Build initialization chains (shared across all regions)
   // -------------------------------------------------------------------------
   // These are built immediately since they don't need temporal operators

   // Only build LatCell binning once (first region)
   static bool BinningBuilt = false;
   if (!BinningBuilt) {
      std::string BinningChainName = "LatCell_CoordinateBinning";

      // Build binning config using makeOpConfig/opParam helpers
      Config BinningConfig =
          makeOpConfig(opParam("NumBins", NumBins), // Number of latitude bins
                       opParam("MinBin", MinLat), // Minimum latitude (radians)
                       opParam("MaxBin", MaxLat)  // Maximum latitude (radians)
          );

      AnalysisManager->parseChainAndBuildOps(BinningChainName, BinningConfig);
      BinningBuilt = true;
   }

   // Return the complete chain string as the stem for temporal operators
   return ChainStr;

} // end buildMOCChain

//------------------------------------------------------------------------------
// Builds the complete operator chain string for computing transect-based MOC
// streamfunction. The chain consists of:
//   Runtime (each timestep):
//     1. PseudoThickness_PseudoToGeometric: Convert to geometric layer
//     thickness
//     2. LayerThickness_BinaryMultiply(NormalVelocity): thickness × velocity
//     3. EdgeTransport_BinaryMultiply(DvEdge): transport × edge width
//     4. TransportField_TransectAccumulator(TransectName): Accumulate across
//     transect
//     5. TransectTransport_PrefixSum: Vertical integration (bottom to top)
//     6. TransectMOC_ScalarMultiply: Convert to Sverdrups (Sv)
//
// This method constructs the complete chain string that will later be parsed
// by buildTemporalChains(), which appends temporal operators and calls
// parseChainAndBuildOps().
//
// Returns the complete chain string (stem) for temporal operator appending
std::string MOC::buildTransectMOCChain(const std::string &TransectName,
                                       Analysis *AnalysisManager) {

   // -------------------------------------------------------------------------
   // Build the complete operator chain string
   // -------------------------------------------------------------------------
   // The complete chain encodes all transformation steps:
   //
   // PseudoThickness_PseudoToGeometric_BinaryMultiply(NormalVelocity)_
   // BinaryMultiply(DvEdge)_TransectAccumulator(TransectName)_PrefixSum_
   // ScalarMultiply(1e-6)
   //
   // This will be parsed by Analysis::parseChainAndBuildOps() which:
   // 1. Splits by underscores to extract field names and operator types
   // 2. Creates operator instances with appropriate configurations
   // 3. Chains them together in sequence

   // Start with PseudoThickness conversion to geometric layer thickness
   std::string ChainStr = "PseudoThickness_PseudoToGeometric";

   // Multiply by NormalVelocity (layer thickness × velocity)
   ChainStr += "_BinaryMultiply(NormalVelocity)";

   // Multiply by DvEdge (transport × edge width)
   ChainStr += "_BinaryMultiply(DvEdge)";

   // Accumulate across transect using mask and sign fields
   ChainStr += "_TransectAccumulator(" + TransectName + ")";

   // Vertical integration (bottom to top) - note: dimension 0 for 1D output
   ChainStr += "_PrefixSum";

   // Convert to Sverdrups
   ChainStr += "_ScalarMultiply(1.0e-6)";

   // -------------------------------------------------------------------------
   // Store configuration for this chain
   // -------------------------------------------------------------------------
   // Each operator in the chain needs configuration. Build a merged config
   // containing parameters for all operators using makeOpConfig/opParam
   // helpers.

   // Build config with all operator parameters
   Config ChainConfig = makeOpConfig(
       opParam("TransectName",
               TransectName), // TransectAccumulator: transect name
       opParam("Dimension",
               0), // PrefixSum: integrate along vertical dimension
                   // (dim 0 of a 1D array from TransectAccumulatorOp)
       opParam("Reverse", true), // PrefixSum: bottom-to-top integration
       opParam("IOName",
               std::string("MOC_streamfunction_transect_") + TransectName));

   // Append config to vector (indexed in same order as combined
   // RegionList+TransectList/ChainStems)
   ChainConfigs.push_back(ChainConfig);

   // Return the complete chain string as the stem for temporal operators
   return ChainStr;

} // end buildTransectMOCChain

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
