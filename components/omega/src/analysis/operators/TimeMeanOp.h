#ifndef OMEGA_TIMEMEANOP_H
#define OMEGA_TIMEMEANOP_H

//===-- analysis/operators/TimeMeanOp.h - TimeMeanOp ------------*- C++ -*-===//
//
/// \file
/// \brief Defines the TimeMeanOp operator for computing time-averaged mean
///
/// TimeMeanOp computes the time-averaged mean of a field over a specified
/// period (e.g., "1day", "1month"). Unlike spatial operators that reduce fields
/// to scalars, TimeMeanOp preserves the full spatial structure of the input
/// field while averaging temporally. The operator accumulates input values at
/// each timestep and divides by the number of accumulations when the period
/// alarm rings.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field, supporting any rank and memory layout. The output Field has the same
/// dimensions and layout as the input Field. Accumulation occurs in-place in
/// the output array using element-wise addition.
///
/// The temporal averaging period is specified via the "Period" configuration
/// option and is managed by an Alarm that signals when to finalize the average.
/// When the alarm rings, the accumulated sum is divided by the count, the
/// result is written to the output, and accumulation restarts for the next
/// period. The operator maintains state (accumulation count and period flag)
/// across timesteps.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"

namespace OMEGA {

/// TimeMeanOp computes the time-averaged mean of a field over a specified
/// period. The operator accumulates input values at each timestep, and when
/// the period alarm rings, divides by the accumulation count to produce the
/// mean. Output Field has the same dimensions and layout as input but always
/// uses Real type for accuracy. Maintains state across timesteps for
/// accumulation.
template <typename ArrayT> class TimeMeanOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same rank as input but always Real precision
   using OutputArrayT = typename std::conditional<
       ArrayT::rank == 1, Array1D_t<Real>,
       typename std::conditional<
           ArrayT::rank == 2, Array2D_t<Real>,
           typename std::conditional<ArrayT::rank == 3, Array3D_t<Real>,
                                     Array4D_t<Real>>::type>::type>::type;

   /// Constructs a TimeMeanOp operator. Reads averaging period from config,
   /// creates output Field matching input dimensions and layout, allocates
   /// output data array for accumulation, and registers the output Field in
   /// the Field registry. Initializes accumulation state. The output Field
   /// name is constructed as InputName + "_TimeMean" + Period (e.g.,
   /// "Temperature_TimeMean1day").
   TimeMeanOp(const std::vector<std::string>
                  &UpstreamNames, ///< [in] input field names
              Config Options      ///< [in] operator config
              )
       : AnalysisOperator("TimeMean") {

      // Store input field names
      InputNames = UpstreamNames;

      // Read averaging period from configuration (e.g., "1day", "1month")
      std::string AvgPeriod;
      Options.get("Period", AvgPeriod);

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_TimeMean" + AvgPeriod;
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Retrieve input Field to extract metadata
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Get dimension info (number of dims and dimension names)
      auto NDims = InputField->getNumDims();
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // Fetch input metadata
      std::string InputDescr, InputUnits, InputStdName;
      ScalarT InputValidMin, InputValidMax;
      InputField->getMetadata("Description", InputDescr);
      InputField->getMetadata("Units", InputUnits);
      InputField->getMetadata("StdName", InputStdName);

      // Register output Field with same dimensions as input but Real type
      auto OutputField =
          Field::create(OutputNames[0],
                        "Time average of " + InputDescr,   // Description
                        InputUnits,                        // Units
                        InputStdName,                      // Standard name
                        -std::numeric_limits<Real>::max(), // Min valid value
                        std::numeric_limits<Real>::max(),  // Max valid value
                        NDims,                             // Rank
                        DimNames                           // Dimension names
          );

      // Store array size for parallel iteration
      ArraySize = static_cast<I4>(InputData.size());

      // Allocate output data array matching input layout but with Real type
      OutputData = OutputArrayT(OutputNames[0] + "_out", InputData.layout());

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

      // Initialize accumulation state
      NumAccum    = 0;
      PeriodAlarm = nullptr;
      IsNewPeriod = true;

      // Default to the flat (non-columnar) path until initialize() detects a
      // recognized horizontal mesh dimension.
      IsMeshDimension = false;
      NHorizOwned     = 0;

   } // end constructor

   /// Initializes the operator after all Fields exist. Determines the index
   /// space (cells, edges, or vertices) from the input field's horizontal
   /// dimension name and, for rank >= 2 mesh fields, stores the appropriate
   /// MinLayer/MaxLayer arrays from VertCoord. TimeMeanOp does NOT abort on an
   /// unrecognized index space: it simply falls back to the flat parallel
   /// dispatch path.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      constexpr I4 InputRank = ArrayT::rank;

      auto InputField = Field::get(InputNames[0]);
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // For rank 1 the sole dimension is the horizontal index space; for rank
      // >= 2 the horizontal index space is the 2nd-to-last dimension.
      std::string IndexSpaceName;
      if constexpr (InputRank == 1) {
         IndexSpaceName = DimNames[0];
      } else {
         IndexSpaceName = DimNames[InputRank - 2];
      }

      IsMeshDimension =
          (IndexSpaceName == "NCells" || IndexSpaceName == "NEdges" ||
           IndexSpaceName == "NVertices");

      // Non-mesh output: leave IsMeshDimension false and use the flat path.
      if (!IsMeshDimension) {
         return;
      }

      // Store owned horizontal count and, for rank >= 2, the active-layer
      // bounds appropriate to the horizontal location.
      if (IndexSpaceName == "NCells") {
         NHorizOwned = Mesh->NCellsOwned;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerCell;
            MaxLayer = VCoord->MaxLayerCell;
         }
      } else if (IndexSpaceName == "NEdges") {
         NHorizOwned = Mesh->NEdgesOwned;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerEdgeBot;
            MaxLayer = VCoord->MaxLayerEdgeTop;
         }
      } else { // NVertices
         NHorizOwned = Mesh->NVerticesOwned;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerVertexBot;
            MaxLayer = VCoord->MaxLayerVertexTop;
         }
      }

   } // end initialize

   /// Sets the period alarm that signals when to finalize the time average.
   /// Called by Analysis during initialization to provide the alarm associated
   /// with the averaging period specified in the constructor config.
   void setPeriodAlarm(Alarm *Alarm ///< [in] alarm for averaging period
                       ) override {
      PeriodAlarm = Alarm;
   }

   /// Computes the time-averaged mean by accumulating input values. On the
   /// first call after alarm rings (IsNewPeriod), copies input to output and
   /// sets count to 1. On subsequent calls, adds input to accumulated sum and
   /// increments count. When period alarm rings, divides accumulated sum by
   /// count to finalize mean and resets state for next period. Updates
   /// timestamp and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Accumulate input values into output array
      if (IsNewPeriod) {
         // Start of new averaging period: initialize with first input
         NumAccum = 1;
         accumInit(InputData);
         IsNewPeriod = false;

         // If the period alarm rings on the first sample of a period, finalize
         // immediately (mean == current value) and start a new period next
         // call.
         if (PeriodAlarm != nullptr && PeriodAlarm->isRinging()) {
            IsNewPeriod = true;
         }
      } else {
         // Continue accumulation: add input to running sum
         accumAdd(InputData);
         ++NumAccum;

         // Check if period alarm is ringing (time to finalize average)
         bool ShouldFinalize =
             (PeriodAlarm != nullptr && PeriodAlarm->isRinging());

         if (ShouldFinalize) {
            // Finalize: divide accumulated sum by count to get mean
            finalize();

            // Reset state for next averaging period
            IsNewPeriod = true;
         }
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

   /// Initializes the output accumulator with the first input sample of a new
   /// averaging period. For column-structured mesh fields the write is
   /// restricted to active layers [MinLayer, MaxLayer] using hierarchical
   /// parallelism. Otherwise a flat parallel loop over the entire array is
   /// used.
   void accumInit(const ArrayT &InputData) {
      OMEGA_SCOPE(LocOutputData, OutputData);
      constexpr I4 InputRank = ArrayT::rank;

      if (IsMeshDimension && InputRank >= 2) {
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocInputData, InputData);

         if constexpr (InputRank == 2) {
            parallelForOuter(
                "TimeMeanInit2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(IHoriz, K) =
                              static_cast<Real>(LocInputData(IHoriz, K));
                       });
                });
         } else if constexpr (InputRank == 3) {
            const I4 Dim0 = InputData.extent(0);
            parallelForOuter(
                "TimeMeanInit3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(I0, IHoriz, K) =
                              static_cast<Real>(LocInputData(I0, IHoriz, K));
                       });
                });
         }
      } else {
         // Flat fallback for scalars, rank-1, and non-mesh outputs
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] =
                    static_cast<Real>(InputData.data()[FlatIdx]);
             });
      }
   } // end accumInit

   /// Adds an input sample to the running accumulator. Uses the same active
   /// layer bounds and dispatch logic as accumInit so that only the layers
   /// initialized at the start of the period are accumulated.
   void accumAdd(const ArrayT &InputData) {
      OMEGA_SCOPE(LocOutputData, OutputData);
      constexpr I4 InputRank = ArrayT::rank;

      if (IsMeshDimension && InputRank >= 2) {
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocInputData, InputData);

         if constexpr (InputRank == 2) {
            parallelForOuter(
                "TimeMeanAdd2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(IHoriz, K) +=
                              static_cast<Real>(LocInputData(IHoriz, K));
                       });
                });
         } else if constexpr (InputRank == 3) {
            const I4 Dim0 = InputData.extent(0);
            parallelForOuter(
                "TimeMeanAdd3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(I0, IHoriz, K) +=
                              static_cast<Real>(LocInputData(I0, IHoriz, K));
                       });
                });
         }
      } else {
         // Flat fallback for scalars, rank-1, and non-mesh outputs
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] +=
                    static_cast<Real>(InputData.data()[FlatIdx]);
             });
      }
   } // end accumAdd

   /// Divides the accumulated sum by the sample count to form the mean. Uses
   /// the same active layer bounds and dispatch logic as accumInit/accumAdd so
   /// fill values in inactive layers are preserved through finalization.
   void finalize() {
      OMEGA_SCOPE(LocOutputData, OutputData);
      const Real InvNumAccum = 1.0 / static_cast<Real>(NumAccum);
      constexpr I4 InputRank = ArrayT::rank;

      if (IsMeshDimension && InputRank >= 2) {
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);

         if constexpr (InputRank == 2) {
            parallelForOuter(
                "TimeMeanFinal2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(IHoriz, K) *= InvNumAccum;
                       });
                });
         } else if constexpr (InputRank == 3) {
            const I4 Dim0 = LocOutputData.extent(0);
            parallelForOuter(
                "TimeMeanFinal3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(I0, IHoriz, K) *= InvNumAccum;
                       });
                });
         }
      } else {
         // Flat fallback for scalars, rank-1, and non-mesh outputs
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] *= InvNumAccum;
             });
      }
   } // end finalize

 private:
   /// Output data array matching input layout but with Real type, used to
   /// accumulate sum and store the final time-averaged mean
   OutputArrayT OutputData;

   /// Number of values accumulated in current averaging period
   I4 NumAccum;

   /// Total size of input/output arrays for the flat parallel path
   I4 ArraySize;

   /// Number of owned points in the horizontal dimension (for MPI-correct
   /// hierarchical loops over mesh-distributed fields)
   I4 NHorizOwned;

   /// Whether the horizontal dimension is a recognized mesh dimension
   /// (NCells/NEdges/NVertices). When false, the flat fallback path is used.
   bool IsMeshDimension;

   /// Min active layer index for each horizontal point (rank >= 2 mesh fields)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (rank >= 2 mesh fields)
   Array1DI4 MaxLayer;

   /// Alarm that signals when averaging period is complete (e.g., end of
   /// day/month)
   Alarm *PeriodAlarm;

   /// Flag indicating whether the next compute() should start a new averaging
   /// period
   bool IsNewPeriod;

}; // end class TimeMeanOp

} // end namespace OMEGA

#endif
