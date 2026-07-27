//===-- Test 5.1: Multi-type operator correctness tests ---------*- C++ -*-===//
//
// Comprehensive unit tests for Analysis operators across supported array types
//
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisOpFactory.h"
#include "Decomp.h"
#include "Field.h"
#include "Forcing.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "TimeStepper.h"
#include "VertCoord.h"
#include "VertMix.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <type_traits>
#include <vector>

using namespace OMEGA;

// Test result tracking
int NumTests  = 0;
int NumPassed = 0;
int NumFailed = 0;

//===----------------------------------------------------------------------===//
// Generic Helper Template Struct
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Template struct consolidating all test helper functions
template <typename ArrayType> struct TestHelper {
   using ScalarT             = typename ArrayType::non_const_value_type;
   static constexpr int Rank = ArrayType::rank;

   // Type-aware tolerance for floating point comparisons
   static ScalarT getTolerance() {
      if constexpr (std::is_integral_v<ScalarT>) {
         return 0; // Exact equality for integers
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         return 1.0e-4f; // Single precision tolerance
      } else {
         return 1.0e-8; // Double precision tolerance
      }
   }

   // Get dimensions based on rank
   static std::vector<I4> getDims(const HorzMesh *Mesh,
                                  const VertCoord *VCoord) {
      if constexpr (Rank == 1) {
         return {Mesh->NCellsSize}; // 1D horizontal array over cells
      } else if constexpr (Rank == 2) {
         return {Mesh->NCellsSize, VCoord->NVertLayers};
      } else if constexpr (Rank == 3) {
         return {Tracers::getNumTracers(), Mesh->NCellsSize,
                 VCoord->NVertLayers};
      }
      return {};
   }

   // Get dimension names
   static std::vector<std::string> getDimNames() {
      if constexpr (Rank == 1) {
         return {"NCells"};
      } else if constexpr (Rank == 2) {
         return {"NCells", "NVertLayers"};
      } else if constexpr (Rank == 3) {
         return {"NTracers", "NCells", "NVertLayers"};
      }
      return {};
   }

   // Create test field for 1D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 1, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4)> ValueFunc) {

      auto DimNames = getDimNames();

      // Set ValidMin/ValidMax based on scalar type to avoid bad_any_cast
      ScalarT ValidMin, ValidMax;
      if constexpr (std::is_integral_v<ScalarT>) {
         ValidMin = std::numeric_limits<ScalarT>::lowest();
         ValidMax = std::numeric_limits<ScalarT>::max();
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         ValidMin = -1.0e30f;
         ValidMax = 1.0e30f;
      } else {
         ValidMin = -1.0e30;
         ValidMax = 1.0e30;
      }

      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", ValidMin, ValidMax, 1, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0]);
      TestField->template attachData<ArrayType>(TestData, false);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         TestDataHost(i) = ValueFunc(i);
      }
      deepCopy(TestData, TestDataHost);
   }

   // Create test field for 2D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 2, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();

      // Set ValidMin/ValidMax based on scalar type to avoid bad_any_cast
      ScalarT ValidMin, ValidMax;
      if constexpr (std::is_integral_v<ScalarT>) {
         ValidMin = std::numeric_limits<ScalarT>::lowest();
         ValidMax = std::numeric_limits<ScalarT>::max();
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         ValidMin = -1.0e30f;
         ValidMax = 1.0e30f;
      } else {
         ValidMin = -1.0e30;
         ValidMax = 1.0e30;
      }

      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", ValidMin, ValidMax, 2, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
      TestField->template attachData<ArrayType>(TestData, false);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            TestDataHost(i, j) = ValueFunc(i, j);
         }
      }
      deepCopy(TestData, TestDataHost);
   }

   // Create test field for 3D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 3, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();

      // Set ValidMin/ValidMax based on scalar type to avoid bad_any_cast
      ScalarT ValidMin, ValidMax;
      if constexpr (std::is_integral_v<ScalarT>) {
         ValidMin = std::numeric_limits<ScalarT>::lowest();
         ValidMax = std::numeric_limits<ScalarT>::max();
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         ValidMin = -1.0e30f;
         ValidMax = 1.0e30f;
      } else {
         ValidMin = -1.0e30;
         ValidMax = 1.0e30;
      }

      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", ValidMin, ValidMax, 3, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1], Dims[2]);
      TestField->template attachData<ArrayType>(TestData, false);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            for (I4 k = 0; k < Dims[2]; ++k) {
               TestDataHost(i, j, k) = ValueFunc(i, j, k);
            }
         }
      }
      deepCopy(TestData, TestDataHost);
   }
};

//------------------------------------------------------------------------------
// Helper function to report test results
void reportTest(const std::string &TestName, bool Passed) {
   NumTests++;
   if (Passed) {
      NumPassed++;
   } else {
      NumFailed++;
      LOG_ERROR("FAIL: {}", TestName);
   }
}

//===----------------------------------------------------------------------===//
// Operator Test Templates
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Template for testing SpatialMaxOp with any array type
template <typename ArrayType>
void testSpatialMaxOpType(const std::string &TypeName, const MachEnv *Env,
                          const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldMax_" + TypeName;

   // Create test field with values based on global cell IDs for MPI correctness
   // Get global cell IDs to ensure unique values across all ranks
   auto Decomp  = Decomp::getDefault();
   auto CellIDH = Decomp->CellIDH; // Global cell IDs (1-based)

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [CellIDH](I4 i) -> ScalarT {
         return static_cast<ScalarT>(CellIDH(i) - 1); // Convert to 0-based
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims,
                          [CellIDH, VCoord](I4 i, I4 j) -> ScalarT {
                             // Unique value: cellID * NVertLayers + layerIndex
                             return static_cast<ScalarT>(
                                 (CellIDH(i) - 1) * VCoord->NVertLayers + j);
                          });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [CellIDH, VCoord, Mesh](I4 i, I4 j, I4 k) -> ScalarT {
             // Unique value: tracerIdx * (NCellsGlobal * NVertLayers) + cellID
             // * NVertLayers + layerIdx
             return static_cast<ScalarT>(
                 i * (Mesh->NCellsGlobal * VCoord->NVertLayers) +
                 (CellIDH(j) - 1) * VCoord->NVertLayers + k);
          });
   }

   // Compute expected max. The operator performs a global reduction across all
   // ranks via MPI_Allreduce, so the expected maximum is based on global mesh
   // properties. The maximum value corresponds to the largest indices in each
   // dimension.
   ScalarT ExpectedMax = 0;

   if constexpr (Rank == 1) {
      // 1D: max is just the highest cell ID (0-based)
      ExpectedMax = static_cast<ScalarT>(Mesh->NCellsGlobal - 1);
   } else if constexpr (Rank == 2) {
      // 2D: max = (NCellsGlobal - 1) * NVertLayers + (NVertLayers - 1)
      ExpectedMax =
          static_cast<ScalarT>((Mesh->NCellsGlobal - 1) * VCoord->NVertLayers +
                               (VCoord->NVertLayers - 1));
   } else if constexpr (Rank == 3) {
      // 3D: max = (NTracers - 1) * (NCellsGlobal * NVertLayers) +
      //           (NCellsGlobal - 1) * NVertLayers + (NVertLayers - 1)
      I4 NTracers = Tracers::getNumTracers();
      ExpectedMax = static_cast<ScalarT>(
          (NTracers - 1) * (Mesh->NCellsGlobal * VCoord->NVertLayers) +
          (Mesh->NCellsGlobal - 1) * VCoord->NVertLayers +
          (VCoord->NVertLayers - 1));
   }

   // Create and compute operator
   auto MaxOp =
       AnalysisOpFactory::createOp("SpatialMax", {FieldName}, makeOpConfig());
   MaxOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MaxOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMax");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   Real ComputedMax     = static_cast<Real>(ResultHost(0));
   Real ExpectedMaxReal = static_cast<Real>(ExpectedMax);

   // Verify
   bool Passed = (std::abs(ComputedMax - ExpectedMaxReal) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMaxOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedMaxReal, ComputedMax);
   }
}

//------------------------------------------------------------------------------
// Template for testing SpatialMinOp with any array type
template <typename ArrayType>
void testSpatialMinOpType(const std::string &TypeName, const MachEnv *Env,
                          const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldMin_" + TypeName;

   // Create test field with values based on global cell IDs for MPI correctness
   auto Decomp  = Decomp::getDefault();
   auto CellIDH = Decomp->CellIDH; // Global cell IDs (1-based)

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [CellIDH](I4 i) -> ScalarT {
         return static_cast<ScalarT>(CellIDH(i) - 1 +
                                     100); // Convert to 0-based, offset by 100
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(
          FieldName, Dims, [CellIDH, VCoord](I4 i, I4 j) -> ScalarT {
             // Unique value: cellID * NVertLayers + layerIndex + offset
             return static_cast<ScalarT>(
                 (CellIDH(i) - 1) * VCoord->NVertLayers + j + 100);
          });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [CellIDH, VCoord, Mesh](I4 i, I4 j, I4 k) -> ScalarT {
             // Unique value: tracerIdx * (NCellsGlobal * NVertLayers) + cellID
             // * NVertLayers + layerIdx + offset
             return static_cast<ScalarT>(
                 i * (Mesh->NCellsGlobal * VCoord->NVertLayers) +
                 (CellIDH(j) - 1) * VCoord->NVertLayers + k + 100);
          });
   }

   // Compute expected min. The operator performs a global reduction across all
   // ranks via MPI_Allreduce, so the expected minimum is based on global mesh
   // properties. The minimum value is always at i=0, j=0 (cell with global ID
   // 0), k=0, plus offset 100.
   ScalarT ExpectedMin = static_cast<ScalarT>(100);

   // Create and compute operator
   auto MinOp =
       AnalysisOpFactory::createOp("SpatialMin", {FieldName}, makeOpConfig());
   MinOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MinOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMin");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   Real ComputedMin     = static_cast<Real>(ResultHost(0));
   Real ExpectedMinReal = static_cast<Real>(ExpectedMin);

   // Verify
   bool Passed = (std::abs(ComputedMin - ExpectedMinReal) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMinOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedMinReal, ComputedMin);
   }
}

//------------------------------------------------------------------------------
// Template for testing SpatialMeanOp with any array type
template <typename ArrayType>
void testSpatialMeanOpType(const std::string &TypeName, const MachEnv *Env,
                           const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldMean_" + TypeName;

   // Create test field with alternating values to properly test mean
   // calculation
   ScalarT Value1 = static_cast<ScalarT>(10);
   ScalarT Value2 = static_cast<ScalarT>(20);

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [Value1, Value2](I4 i) -> ScalarT {
         return ((i % 2) == 0) ? Value1 : Value2;
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j) -> ScalarT {
                             return (((i + j) % 2) == 0) ? Value1 : Value2;
                          });
   } else if constexpr (Rank == 3) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j, I4 k) -> ScalarT {
                             return (((i + j + k) % 2) == 0) ? Value1 : Value2;
                          });
   }

   // Calculate expected mean by counting actual Value1 and Value2 elements
   // in the active region (accounting for masked layers)
   I8 LocalCount1 = 0, LocalCount2 = 0;

   if constexpr (Rank == 1) {
      // 1D: count based on horizontal index pattern
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if ((i % 2) == 0)
            LocalCount1++;
         else
            LocalCount2++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: count based on (i+j) pattern within active vertical layers
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (((i + j) % 2) == 0)
               LocalCount1++;
            else
               LocalCount2++;
         }
      }
   } else if constexpr (Rank == 3) {
      // 3D: count based on (t+i+j) pattern within active layers
      I4 NTracers = Tracers::getNumTracers();
      for (I4 t = 0; t < NTracers; ++t) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
            for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
                 ++j) {
               if (((t + i + j) % 2) == 0)
                  LocalCount1++;
               else
                  LocalCount2++;
            }
         }
      }
   }

   // Global sum across MPI ranks
   I8 Count1         = globalSum(LocalCount1, Env->getComm());
   I8 Count2         = globalSum(LocalCount2, Env->getComm());
   I8 TotalElements  = Count1 + Count2;
   Real ExpectedMean = (static_cast<Real>(Value1) * static_cast<Real>(Count1) +
                        static_cast<Real>(Value2) * static_cast<Real>(Count2)) /
                       static_cast<Real>(TotalElements);

   // Create and compute operator
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, makeOpConfig());
   MeanOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialMean");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   Real ComputedMean = ResultHost(0);

   // Verify
   bool Passed = (std::abs(ComputedMean - ExpectedMean) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMeanOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedMean, ComputedMean);
   }
}

//------------------------------------------------------------------------------
// Template for testing SpatialStdDevOp with any array type
template <typename ArrayType>
void testSpatialStdDevOpType(const std::string &TypeName, const MachEnv *Env,
                             const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldStdDev_" + TypeName;

   // Create test field with alternating values to properly test std dev
   // calculation
   ScalarT Value1 = static_cast<ScalarT>(10);
   ScalarT Value2 = static_cast<ScalarT>(20);

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [Value1, Value2](I4 i) -> ScalarT {
         return ((i % 2) == 0) ? Value1 : Value2;
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j) -> ScalarT {
                             return (((i + j) % 2) == 0) ? Value1 : Value2;
                          });
   } else if constexpr (Rank == 3) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j, I4 k) -> ScalarT {
                             return (((i + j + k) % 2) == 0) ? Value1 : Value2;
                          });
   }

   // Calculate expected standard deviation by counting actual Value1 and Value2
   // elements in the active region (accounting for masked layers)
   I8 LocalCount1 = 0, LocalCount2 = 0;

   if constexpr (Rank == 1) {
      // 1D: count based on horizontal index pattern
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if ((i % 2) == 0)
            LocalCount1++;
         else
            LocalCount2++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: count based on (i+j) pattern within active vertical layers
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (((i + j) % 2) == 0)
               LocalCount1++;
            else
               LocalCount2++;
         }
      }
   } else if constexpr (Rank == 3) {
      // 3D: count based on (t+i+j) pattern within active layers
      I4 NTracers = Tracers::getNumTracers();
      for (I4 t = 0; t < NTracers; ++t) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
            for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
                 ++j) {
               if (((t + i + j) % 2) == 0)
                  LocalCount1++;
               else
                  LocalCount2++;
            }
         }
      }
   }

   // Global sum across MPI ranks
   I8 Count1        = globalSum(LocalCount1, Env->getComm());
   I8 Count2        = globalSum(LocalCount2, Env->getComm());
   I8 TotalElements = Count1 + Count2;
   Real Mean        = (static_cast<Real>(Value1) * static_cast<Real>(Count1) +
                static_cast<Real>(Value2) * static_cast<Real>(Count2)) /
               static_cast<Real>(TotalElements);

   // Standard deviation: sqrt(sum((x_i - mean)^2) / N)
   Real SumSquaredDiff = static_cast<Real>(Count1) *
                             std::pow(static_cast<Real>(Value1) - Mean, 2.0) +
                         static_cast<Real>(Count2) *
                             std::pow(static_cast<Real>(Value2) - Mean, 2.0);
   Real ExpectedStdDev =
       std::sqrt(SumSquaredDiff / static_cast<Real>(TotalElements));

   // SpatialStdDevOp requires a pre-existing _SpatialMean field for the input.
   // Create and compute a SpatialMeanOp first so that field is registered.
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, makeOpConfig());
   MeanOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Now create and compute the StdDev operator
   auto StdDevOp = AnalysisOpFactory::createOp("SpatialStdDev", {FieldName},
                                               makeOpConfig());
   StdDevOp->initialize(Env, Mesh, VCoord, makeOpConfig());
   StdDevOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialStdDev");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   Real ComputedStdDev = ResultHost(0);

   // Verify
   bool Passed = (std::abs(ComputedStdDev - ExpectedStdDev) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialStdDevOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedStdDev, ComputedStdDev);
   }
}

//------------------------------------------------------------------------------
// Template for testing TimeMeanOp with any array type
template <typename ArrayType>
void testTimeMeanOpType(const std::string &TypeName, const MachEnv *Env,
                        const HorzMesh *Mesh, const VertCoord *VCoord,
                        Clock *ModelClock) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldTimeMean_" + TypeName;

   // Create test field with initial value
   ScalarT BaseValue = static_cast<ScalarT>(5);

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims,
                          [BaseValue](I4 i) -> ScalarT { return BaseValue; });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims, [BaseValue](I4 i, I4 j) -> ScalarT {
         return BaseValue;
      });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [BaseValue](I4 i, I4 j, I4 k) -> ScalarT { return BaseValue; });
   }

   // Get the field and its data array for updating during time loop
   auto TestField = Field::get(FieldName);
   auto TestData  = TestField->template getDataArray<ArrayType>();

   // Set up time stepping parameters
   const int NumSteps        = 5; // Accumulate over 5 timesteps
   TimeInterval StepInterval = ModelClock->getTimeStep();

   // Calculate period interval (NumSteps * timestep)
   R8 StepSeconds;
   StepInterval.get(StepSeconds, TimeUnits::Seconds);
   R8 PeriodSeconds = StepSeconds * NumSteps;
   TimeInterval PeriodInterval(PeriodSeconds, TimeUnits::Seconds);

   // Create TimeMeanOp with a valid period string (e.g., "5seconds")
   // The Period string is just a label used in the output field name
   Config OpConfig;
   std::string PeriodLabel =
       std::to_string(static_cast<int>(PeriodSeconds)) + "seconds";
   //   OpConfig.set("Period", PeriodLabel);

   auto TimeMeanOp = AnalysisOpFactory::createOp(
       "TimeMean", {FieldName}, makeOpConfig(opParam("Period", PeriodLabel)));
   TimeMeanOp->initialize(Env, Mesh, VCoord, OpConfig);

   // Create a period alarm that rings after NumSteps
   TimeInstant StartTime = ModelClock->getCurrentTime();
   Alarm PeriodAlarm("TestPeriodAlarm_" + TypeName, PeriodInterval, StartTime);
   TimeMeanOp->setPeriodAlarm(&PeriodAlarm);

   // Time-stepping loop: update field values and compute mean at each step
   std::vector<ScalarT> ValuesAtEachStep;

   for (int step = 0; step < NumSteps; ++step) {
      // Update field values to simulate time evolution
      // Value at each step = BaseValue + step (e.g., 5, 6, 7, 8, 9)
      ScalarT CurrentValue = static_cast<ScalarT>(static_cast<Real>(BaseValue) +
                                                  static_cast<Real>(step));
      ValuesAtEachStep.push_back(CurrentValue);

      // Update the field data on device
      auto TestDataHost = createHostMirrorCopy(TestData);

      if constexpr (Rank == 1) {
         for (I4 i = 0; i < Dims[0]; ++i) {
            TestDataHost(i) = CurrentValue;
         }
      } else if constexpr (Rank == 2) {
         for (I4 i = 0; i < Dims[0]; ++i) {
            for (I4 j = 0; j < Dims[1]; ++j) {
               TestDataHost(i, j) = CurrentValue;
            }
         }
      } else if constexpr (Rank == 3) {
         for (I4 i = 0; i < Dims[0]; ++i) {
            for (I4 j = 0; j < Dims[1]; ++j) {
               for (I4 k = 0; k < Dims[2]; ++k) {
                  TestDataHost(i, j, k) = CurrentValue;
               }
            }
         }
      }

      deepCopy(TestData, TestDataHost);

      // Advance clock to next timestep
      ModelClock->advance();
      TimeInstant CurrentTime = ModelClock->getCurrentTime();

      // Update alarm status based on current time
      PeriodAlarm.updateStatus(CurrentTime);

      // Compute the time mean (accumulates internally)
      TimeMeanOp->compute(CurrentTime);

      // Check if alarm is ringing (should ring after last step)
      if (PeriodAlarm.isRinging()) {
         // Mean should now be finalized
         break;
      }
   }

   // Calculate expected mean: average of [BaseValue, BaseValue+1, ...,
   // BaseValue+(NumSteps-1)] For BaseValue=5 and NumSteps=5: avg of [5, 6, 7,
   // 8, 9] = 7.0
   Real Sum = 0.0;
   for (const auto &val : ValuesAtEachStep) {
      Sum += static_cast<Real>(val);
   }
   Real ExpectedMean = Sum / static_cast<Real>(NumSteps);

   // Get result field - output field name includes the Period label
   std::string ResultFieldName = FieldName + "_TimeMean" + PeriodLabel;
   auto ResultField            = Field::get(ResultFieldName);

   // Verify a sample of values. The TimeMeanOp output field is always Real type
   // regardless of input type.
   bool Passed = true;
   if constexpr (Rank == 1) {
      auto ResultData = ResultField->getDataArray<Array1D_t<Real>>();
      auto ResultHost = createHostMirrorCopy(ResultData);

      for (I4 i = 0; i < std::min(10, Dims[0]); ++i) {
         Real ComputedValue = ResultHost(i);
         if (std::abs(ComputedValue - ExpectedMean) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            LOG_ERROR("  At index {}: Expected {}, Got {}", i, ExpectedMean,
                      ComputedValue);
            break;
         }
      }
   } else if constexpr (Rank == 2) {
      auto ResultData = ResultField->getDataArray<Array2D_t<Real>>();
      auto ResultHost = createHostMirrorCopy(ResultData);

      for (I4 i = 0; i < std::min(5, Dims[0]); ++i) {
         for (I4 j = 0; j < std::min(5, Dims[1]); ++j) {
            Real ComputedValue = ResultHost(i, j);
            if (std::abs(ComputedValue - ExpectedMean) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               LOG_ERROR("  At index ({}, {}): Expected {}, Got {}", i, j,
                         ExpectedMean, ComputedValue);
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      auto ResultData = ResultField->getDataArray<Array3D_t<Real>>();
      auto ResultHost = createHostMirrorCopy(ResultData);

      for (I4 i = 0; i < std::min(3, Dims[0]); ++i) {
         for (I4 j = 0; j < std::min(3, Dims[1]); ++j) {
            for (I4 k = 0; k < std::min(3, Dims[2]); ++k) {
               Real ComputedValue = ResultHost(i, j, k);
               if (std::abs(ComputedValue - ExpectedMean) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
                  LOG_ERROR("  At index ({}, {}, {}): Expected {}, Got {}", i,
                            j, k, ExpectedMean, ComputedValue);
                  break;
               }
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   reportTest("TimeMeanOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected mean: {}", ExpectedMean);
      LOG_ERROR("  Period label: {}", PeriodLabel);
   }
}

//===----------------------------------------------------------------------===//
// Main Test Functions
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Test SpatialMaxOp with all array types
void testSpatialMaxOp(const MachEnv *Env, const HorzMesh *Mesh,
                      const VertCoord *VCoord) {

   // 1D arrays - 4 scalar types
   testSpatialMaxOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays - 4 scalar types
   testSpatialMaxOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays - 4 scalar types
   testSpatialMaxOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test SpatialMinOp with all array types
void testSpatialMinOp(const MachEnv *Env, const HorzMesh *Mesh,
                      const VertCoord *VCoord) {

   // 1D arrays
   testSpatialMinOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialMinOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testSpatialMinOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialMinOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testSpatialMinOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialMinOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test SpatialMeanOp with all array types
void testSpatialMeanOp(const MachEnv *Env, const HorzMesh *Mesh,
                       const VertCoord *VCoord) {

   // 1D arrays
   testSpatialMeanOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testSpatialMeanOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testSpatialMeanOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test SpatialStdDevOp with all array types
void testSpatialStdDevOp(const MachEnv *Env, const HorzMesh *Mesh,
                         const VertCoord *VCoord) {

   // 1D arrays
   testSpatialStdDevOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testSpatialStdDevOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testSpatialStdDevOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test TimeMeanOp with all array types
void testTimeMeanOp(const MachEnv *Env, const HorzMesh *Mesh,
                    const VertCoord *VCoord, Clock *ModelClock) {

   // 1D arrays
   testTimeMeanOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord, ModelClock);

   // 2D arrays
   testTimeMeanOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord, ModelClock);

   // 3D arrays
   testTimeMeanOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord, ModelClock);
}

//------------------------------------------------------------------------------
// Template for testing ScalarMultiplyOp with any array type
template <typename ArrayType>
void testScalarMultiplyOpType(const std::string &TypeName, const MachEnv *Env,
                              const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldScalarMult_" + TypeName;

   // Create test field with known values
   ScalarT BaseValue = static_cast<ScalarT>(10);
   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims,
                          [BaseValue](I4 i) -> ScalarT { return BaseValue; });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims, [BaseValue](I4 i, I4 j) -> ScalarT {
         return BaseValue;
      });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [BaseValue](I4 i, I4 j, I4 k) -> ScalarT { return BaseValue; });
   }
   // Test scalar multiplication
   std::string ScalarStr = "2.5";
   Real Scalar           = 2.5;
   ScalarT ExpectedValue =
       static_cast<ScalarT>(static_cast<Real>(BaseValue) * Scalar);
   // Create and compute operator
   Config OpConfig;
   OpConfig.add("Scalar", Scalar);

   auto MultOp =
       AnalysisOpFactory::createOp("ScalarMultiply", {FieldName}, OpConfig);
   MultOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MultOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(FieldName + "_ScalarMultiply(" + ScalarStr + ")");
   auto ResultData = ResultField->getDataArray<ArrayType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify all owned cells
   bool Passed = true;
   if constexpr (Rank == 1) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if (std::abs(static_cast<Real>(ResultHost(i)) -
                      static_cast<Real>(ExpectedValue)) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            break;
         }
      }
   } else if constexpr (Rank == 2) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (std::abs(static_cast<Real>(ResultHost(i, j)) -
                         static_cast<Real>(ExpectedValue)) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            for (I4 k = VCoord->MinLayerCellH(j); k <= VCoord->MaxLayerCellH(j);
                 ++k) {
               if (std::abs(static_cast<Real>(ResultHost(i, j, k)) -
                            static_cast<Real>(ExpectedValue)) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
                  break;
               }
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   reportTest("ScalarMultiplyOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Scalar: {}", static_cast<Real>(ExpectedValue),
                Scalar);
   }
}

//------------------------------------------------------------------------------
// Test ScalarMultiplyOp with all array types
void testScalarMultiplyOp(const MachEnv *Env, const HorzMesh *Mesh,
                          const VertCoord *VCoord) {

   // 1D arrays
   testScalarMultiplyOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testScalarMultiplyOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testScalarMultiplyOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Template for testing BinaryMultiplyOp with same-rank arrays
template <typename ArrayType>
void testBinaryMultiplyOpSameRank(const std::string &TypeName,
                                  const MachEnv *Env, const HorzMesh *Mesh,
                                  const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims   = Helper::getDims(Mesh, VCoord);
   std::string Field1Name = "TestField1BinaryMult_" + TypeName;
   std::string Field2Name = "TestField2BinaryMult_" + TypeName;

   // Create two test fields with known values
   ScalarT Value1 = static_cast<ScalarT>(3);
   ScalarT Value2 = static_cast<ScalarT>(4);

   if constexpr (Rank == 1) {
      Helper::createField(Field1Name, Dims,
                          [Value1](I4 i) -> ScalarT { return Value1; });
      Helper::createField(Field2Name, Dims,
                          [Value2](I4 i) -> ScalarT { return Value2; });
   } else if constexpr (Rank == 2) {
      Helper::createField(Field1Name, Dims,
                          [Value1](I4 i, I4 j) -> ScalarT { return Value1; });
      Helper::createField(Field2Name, Dims,
                          [Value2](I4 i, I4 j) -> ScalarT { return Value2; });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          Field1Name, Dims,
          [Value1](I4 i, I4 j, I4 k) -> ScalarT { return Value1; });
      Helper::createField(
          Field2Name, Dims,
          [Value2](I4 i, I4 j, I4 k) -> ScalarT { return Value2; });
   }

   ScalarT ExpectedValue = static_cast<ScalarT>(static_cast<Real>(Value1) *
                                                static_cast<Real>(Value2));

   // Create and compute operator
   auto MultOp = AnalysisOpFactory::createOp(
       "BinaryMultiply", {Field1Name, Field2Name}, makeOpConfig());
   MultOp->initialize(Env, Mesh, VCoord, makeOpConfig());
   TimeInstant TestTime;
   MultOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(Field1Name + "_BinaryMultiply(" + Field2Name + ")");
   auto ResultData = ResultField->getDataArray<ArrayType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify all owned cells
   bool Passed = true;
   if constexpr (Rank == 1) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if (std::abs(static_cast<Real>(ResultHost(i)) -
                      static_cast<Real>(ExpectedValue)) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            break;
         }
      }
   } else if constexpr (Rank == 2) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (std::abs(static_cast<Real>(ResultHost(i, j)) -
                         static_cast<Real>(ExpectedValue)) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            for (I4 k = VCoord->MinLayerCellH(j); k <= VCoord->MaxLayerCellH(j);
                 ++k) {
               if (std::abs(static_cast<Real>(ResultHost(i, j, k)) -
                            static_cast<Real>(ExpectedValue)) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
                  break;
               }
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   reportTest("BinaryMultiplyOp (same-rank): " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Value1: {}, Value2: {}",
                static_cast<Real>(ExpectedValue), static_cast<Real>(Value1),
                static_cast<Real>(Value2));
   }
}

//------------------------------------------------------------------------------
// Test BinaryMultiplyOp with vertical expansion (2D x 1D)
template <typename Array2DType>
void testBinaryMultiplyOpVerticalExpansion(const std::string &TypeName,
                                           const MachEnv *Env,
                                           const HorzMesh *Mesh,
                                           const VertCoord *VCoord) {

   using Helper2D       = TestHelper<Array2DType>;
   using ScalarT        = typename Helper2D::ScalarT;
   constexpr int Rank2D = Helper2D::Rank;

   // Only test 2D arrays for vertical expansion
   if constexpr (Rank2D != 2) {
      return;
   }

   std::vector<I4> Dims2D = Helper2D::getDims(Mesh, VCoord);
   std::vector<I4> Dims1D = {Dims2D[0]}; // Just horizontal dimension

   std::string Field2DName = "TestField2DVertExp_" + TypeName;
   std::string Field1DName = "TestField1DVertExp_" + TypeName;

   // Create 2D field with known values
   ScalarT Value2D = static_cast<ScalarT>(5);
   Helper2D::createField(Field2DName, Dims2D,
                         [Value2D](I4 i, I4 j) -> ScalarT { return Value2D; });

   // Create 1D field (horizontal only)
   ScalarT Value1D = static_cast<ScalarT>(2);
   using Helper1D  = TestHelper<Array1DR8>; // Use same scalar type base
   Helper1D::createField(Field1DName, Dims1D,
                         [Value1D](I4 i) -> ScalarT { return Value1D; });

   ScalarT ExpectedValue = static_cast<ScalarT>(static_cast<Real>(Value2D) *
                                                static_cast<Real>(Value1D));

   // Create and compute operator (2D x 1D = 2D with vertical expansion)
   auto MultOp = AnalysisOpFactory::createOp(
       "BinaryMultiply", {Field2DName, Field1DName}, makeOpConfig());
   MultOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MultOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(Field2DName + "_BinaryMultiply(" + Field1DName + ")");
   auto ResultData = ResultField->getDataArray<Array2DType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify: all active vertical layers should have same value (1D value
   // replicated)
   bool Passed = true;
   for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
      for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
           ++j) {
         if (std::abs(static_cast<Real>(ResultHost(i, j)) -
                      static_cast<Real>(ExpectedValue)) >
             static_cast<Real>(Helper2D::getTolerance())) {
            Passed = false;
            LOG_ERROR("  At ({}, {}): Expected {}, Got {}", i, j,
                      static_cast<Real>(ExpectedValue),
                      static_cast<Real>(ResultHost(i, j)));
            break;
         }
      }
      if (!Passed)
         break;
   }

   reportTest("BinaryMultiplyOp (vertical expansion): " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Value2D: {}, Value1D: {}",
                static_cast<Real>(ExpectedValue), static_cast<Real>(Value2D),
                static_cast<Real>(Value1D));
   }
}

//------------------------------------------------------------------------------
// Test BinaryMultiplyOp with all array types
void testBinaryMultiplyOp(const MachEnv *Env, const HorzMesh *Mesh,
                          const VertCoord *VCoord) {

   // Test same-rank multiplication for 1D arrays
   testBinaryMultiplyOpSameRank<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // Test same-rank multiplication for 2D arrays
   testBinaryMultiplyOpSameRank<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // Test same-rank multiplication for 3D arrays
   testBinaryMultiplyOpSameRank<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array3DR8>("3D-R8", Env, Mesh, VCoord);

   // Test vertical expansion for 2D × 1D
   testBinaryMultiplyOpVerticalExpansion<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpVerticalExpansion<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpVerticalExpansion<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpVerticalExpansion<Array2DR8>("2D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Template for testing PrefixSumOp with 2D arrays
template <typename ArrayType>
void testPrefixSumOpType(const std::string &TypeName, const MachEnv *Env,
                         const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;
   // Only test 2D arrays for PrefixSum (most relevant for MOC)
   if constexpr (Rank != 2) {
      return;
   }

   // Use bins (non-distributed dimension) for testing
   // Create a simple 2D field: NBins x NVertLayers
   const I4 NBins       = 10; // Small number for testing
   const I4 NVertLayers = VCoord->NVertLayers;

   std::string FieldName = "TestFieldPrefixSum_" + TypeName;

   Dimension::create("TestBins", NBins);
   // Create field with bins × vertical levels
   auto TestField =
       Field::create(FieldName, "Test field for PrefixSum", "units", "",
                     static_cast<ScalarT>(-1000), static_cast<ScalarT>(1000), 2,
                     {"TestBins", "NVertLayers"});
   // Allocate and initialize data: value = (bin_idx + 1) * (level + 1)
   ArrayType TestData(FieldName + "_data", NBins, NVertLayers);
   auto TestHost = Kokkos::create_mirror_view(TestData);
   for (I4 i = 0; i < NBins; ++i) {
      for (I4 j = 0; j < NVertLayers; ++j) {
         TestHost(i, j) = static_cast<ScalarT>((i + 1) * (j + 1));
      }
   }
   deepCopy(TestData, TestHost);
   TestField->template attachData<ArrayType>(TestData, false);
   // Test vertical integration (dimension 1, reverse = true for MOC)
   auto PrefixOp = AnalysisOpFactory::createOp(
       "PrefixSum", {FieldName},
       makeOpConfig(opParam("Dimension", 1), opParam("Reverse", true)));
   PrefixOp->initialize(Env, Mesh, VCoord, makeOpConfig());
   TimeInstant TestTime;
   PrefixOp->compute(TestTime);
   // Get result
   auto ResultField = Field::get(FieldName + "_PrefixSum");
   auto ResultData  = ResultField->getDataArray<ArrayType>();
   auto ResultHost  = createHostMirrorCopy(ResultData);
   // Verify reverse cumulative sum
   bool Passed = true;

   for (I4 i = 0; i < NBins; ++i) {
      for (I4 j = 0; j < NVertLayers; ++j) {
         // Compute expected: sum from j to NVertLayers-1
         Real Expected = 0.0;
         for (I4 k = j; k < NVertLayers; ++k) {
            Expected += static_cast<Real>((i + 1) * (k + 1));
         }

         Real Computed = static_cast<Real>(ResultHost(i, j));

         if (std::abs(Computed - Expected) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            LOG_ERROR("  At bin {} level {}: Expected {}, Got {}", i, j,
                      Expected, Computed);
            break;
         }
      }
      if (!Passed)
         break;
   }
   reportTest("PrefixSumOp: " + TypeName, Passed);
}

//------------------------------------------------------------------------------
// Test PrefixSumOp with 2D array types
void testPrefixSumOp(const MachEnv *Env, const HorzMesh *Mesh,
                     const VertCoord *VCoord) {

   // Test 2D arrays (most relevant for MOC vertical integration)
   testPrefixSumOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testPrefixSumOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testPrefixSumOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testPrefixSumOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test CoordinateBinningOp with 1D coordinate arrays
void testCoordinateBinningOp(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   std::string FieldName = "TestCoordBinning";

   // Create a simple 1D coordinate field with values from -90 to 90 (like
   // latitude)
   I4 NCells = Mesh->NCellsSize;
   auto TestField =
       Field::create(FieldName, "Test coordinate field for binning", "degrees",
                     "", -180.0, 180.0, 1, {"NCells"});

   Array1DR8 TestData(FieldName + "_data", NCells);
   TestField->attachData<Array1DR8>(TestData, false);

   auto TestDataHost = Kokkos::create_mirror_view(TestData);
   auto Decomp       = Decomp::getDefault();
   auto CellIDH      = Decomp->CellIDH;

   // Create latitude-like distribution: map global cell ID to latitude
   for (I4 i = 0; i < NCells; ++i) {
      Real Fraction = static_cast<Real>(CellIDH(i) - 1) /
                      static_cast<Real>(Mesh->NCellsGlobal - 1);
      TestDataHost(i) = -90.0 + Fraction * 180.0; // -90 to +90
   }
   deepCopy(TestData, TestDataHost);

   // Create CoordinateBinning operator
   I4 NumBins     = 18; // 10-degree bins
   auto BinningOp = AnalysisOpFactory::createOp(
       "CoordinateBinning", {FieldName},
       makeOpConfig(opParam("NumBins", NumBins), opParam("MinBin", -90.0),
                    opParam("MaxBin", 90.0)));
   BinningOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   BinningOp->compute(TestTime);

   // Get result (bin indices)
   auto ResultField = Field::get(FieldName + "_BinIndex");
   auto ResultData  = ResultField->getDataArray<Array1DI4>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Verify bin assignments
   bool Passed = true;
   for (I4 i = 0; i < std::min(10, NCells); ++i) {
      Real Coord     = TestDataHost(i);
      I4 ComputedBin = ResultHost(i);
      Real BinWidth  = 180.0 / static_cast<Real>(NumBins);
      I4 ExpectedBin = static_cast<I4>((Coord + 90.0) / BinWidth);

      // Handle edge case at max boundary
      if (ExpectedBin >= NumBins)
         ExpectedBin = NumBins - 1;

      if (ComputedBin != ExpectedBin) {
         Passed = false;
         LOG_ERROR("  At index {}: Coord={}, Expected bin {}, Got bin {}", i,
                   Coord, ExpectedBin, ComputedBin);
         break;
      }
   }

   reportTest("CoordinateBinningOp", Passed);
}

//------------------------------------------------------------------------------
// Test BinnedAccumulatorOp with 2D field and bin indices
void testBinnedAccumulatorOp(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   // Create bin index field (assign cells to bins based on pattern)
   I4 NCells  = Mesh->NCellsSize;
   I4 NumBins = 5;

   std::string BinIndexFieldName = "TestBinIndices";
   auto BinField = Field::create(BinIndexFieldName, "Test bin indices", "1", "",
                                 0, NumBins, 1, {"NCells"});

   Array1DI4 BinData(BinIndexFieldName + "_data", NCells);
   BinField->attachData<Array1DI4>(BinData, false);

   auto BinDataHost = Kokkos::create_mirror_view(BinData);
   for (I4 i = 0; i < NCells; ++i) {
      BinDataHost(i) = i % NumBins; // Distribute cells across bins
   }
   deepCopy(BinData, BinDataHost);

   // Create value field to accumulate (2D: NCells × NVertLayers)
   std::string ValueFieldName = "TestAccumValues";
   I4 NVertLayers             = VCoord->NVertLayers;

   auto ValueField =
       Field::create(ValueFieldName, "Test values for accumulation", "m3/s", "",
                     -1e30, 1e30, 2, {"NCells", "NVertLayers"});

   Array2DR8 ValueData(ValueFieldName + "_data", NCells, NVertLayers);
   ValueField->attachData<Array2DR8>(ValueData, false);

   auto ValueDataHost = Kokkos::create_mirror_view(ValueData);
   for (I4 i = 0; i < NCells; ++i) {
      for (I4 j = 0; j < NVertLayers; ++j) {
         ValueDataHost(i, j) = 1.0; // Uniform value for easy verification
      }
   }
   deepCopy(ValueData, ValueDataHost);

   // Create BinnedAccumulator operator
   auto AccumOp = AnalysisOpFactory::createOp(
       "BinnedAccumulator", {ValueFieldName, BinIndexFieldName},
       makeOpConfig(opParam("NumBins", NumBins)));
   AccumOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   AccumOp->compute(TestTime);

   // Get result (should be NumBins × NVertLayers)
   auto ResultField = Field::get(ValueFieldName + "_BinnedAccumulator(" +
                                 BinIndexFieldName + ")");
   auto ResultData  = ResultField->getDataArray<Array2DR8>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Calculate expected sums per bin per layer, respecting MinLayer/MaxLayer
   // (only cells where layer K is active contribute to that layer's bin sum)
   // Local counts per bin per layer
   std::vector<std::vector<I8>> LocalCounts(NumBins,
                                            std::vector<I8>(NVertLayers, 0));
   for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
      I4 BinIdx = BinDataHost(i);
      I4 KMin   = VCoord->MinLayerCellH(i);
      I4 KMax   = VCoord->MaxLayerCellH(i);
      for (I4 k = KMin; k <= KMax; ++k) {
         LocalCounts[BinIdx][k]++;
      }
   }

   // Global sum across MPI ranks
   std::vector<std::vector<I8>> GlobalCounts(NumBins,
                                             std::vector<I8>(NVertLayers, 0));
   for (I4 bin = 0; bin < NumBins; ++bin) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         GlobalCounts[bin][k] = globalSum(LocalCounts[bin][k], Env->getComm());
      }
   }

   // Verify: each bin/layer sum should equal count * 1.0 (uniform value field)
   bool Passed = true;
   for (I4 bin = 0; bin < NumBins; ++bin) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         Real Expected = static_cast<Real>(GlobalCounts[bin][k]) * 1.0;
         Real Computed = ResultHost(bin, k);

         if (std::abs(Computed - Expected) > 1.0e-8) {
            Passed = false;
            LOG_ERROR("  Bin {}, Layer {}: Expected {}, Got {}", bin, k,
                      Expected, Computed);
            break;
         }
      }
      if (!Passed)
         break;
   }

   reportTest("BinnedAccumulatorOp", Passed);
}

//------------------------------------------------------------------------------
// Template for testing PseudoToGeometricOp with any array type
template <typename ArrayType>
void testPseudoToGeometricOpType(const std::string &TypeName,
                                 const MachEnv *Env, const HorzMesh *Mesh,
                                 const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   // Only test Real types (R4, R8) since SpecVol multiplication requires reals
   if constexpr (!std::is_floating_point_v<ScalarT>) {
      return;
   }

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldPseudoToGeom_" + TypeName;

   // Create and populate SpecVol Field with test data
   auto SpecVolField = Field::get("SpecVol");
   if (!SpecVolField) {
      // Create SpecVol field if it doesn't exist
      SpecVolField = Field::create("SpecVol", "Test specific volume", "m^3/kg",
                                   "", 0.0, 10.0, 2, {"NCells", "NVertLayers"});

      Array2DReal SpecVolDataNew("SpecVol_data", Mesh->NCellsSize,
                                 VCoord->NVertLayers);
      SpecVolField->attachData<Array2DReal>(SpecVolDataNew, false);
   }

   auto SpecVolData = SpecVolField->getDataArray<Array2DReal>();
   auto SpecVolHost = Kokkos::create_mirror_view(SpecVolData);

   // Fill with test values: SpecVol(i,k) = 1.0 + 0.01*i + 0.001*k
   for (I4 i = 0; i < Mesh->NCellsSize; ++i) {
      for (I4 k = 0; k < VCoord->NVertLayers; ++k) {
         SpecVolHost(i, k) = 1.0 + 0.01 * i + 0.001 * k;
      }
   }
   deepCopy(SpecVolData, SpecVolHost);

   // Create pseudo-coordinate test field with known values
   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [](I4 i) -> ScalarT {
         return static_cast<ScalarT>(i + 1); // Simple linear values
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims, [](I4 i, I4 k) -> ScalarT {
         return static_cast<ScalarT>((i + 1) * (k + 1)); // Product of indices
      });
   } else if constexpr (Rank == 3) {
      Helper::createField(FieldName, Dims, [](I4 i, I4 j, I4 k) -> ScalarT {
         return static_cast<ScalarT>((i + 1) * (j + 1) * (k + 1));
      });
   }

   // Create and compute operator
   auto PseudoToGeomOp = AnalysisOpFactory::createOp(
       "PseudoToGeometric", {FieldName}, makeOpConfig());
   PseudoToGeomOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   PseudoToGeomOp->compute(TestTime);

   // Get result
   auto ResultField = Field::get(FieldName + "_PseudoToGeometric");
   auto ResultData  = ResultField->getDataArray<ArrayType>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Get input field for verification
   auto InputField = Field::get(FieldName);
   auto InputData  = InputField->getDataArray<ArrayType>();
   auto InputHost  = createHostMirrorCopy(InputData);

   // Verify all elements
   bool Passed     = true;
   Real Tolerance  = static_cast<Real>(Helper::getTolerance());
   int NumVerified = 0;

   if constexpr (Rank == 1) {
      // 1D: check all horizontal points
      for (I4 i = 0; i < Mesh->NCellsOwned; i++) {
         // Compute expected using same precision path as operator
         Real PseudoVal  = static_cast<Real>(InputHost(i));
         Real SpecVolVal = SpecVolHost(i, 0);
         ScalarT Expected =
             static_cast<ScalarT>(RhoSw * SpecVolVal * PseudoVal);
         ScalarT Computed = ResultHost(i);

         if (std::abs(static_cast<Real>(Computed) -
                      static_cast<Real>(Expected)) /
                 static_cast<Real>(Expected) >
             Tolerance) {
            Passed = false;
            break;
         }
         NumVerified++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: check all (i,k) points in active layers
      for (I4 i = 0; i < Mesh->NCellsOwned; i++) {
         // Respect active layer bounds for this cell, clamped to array bounds
         I4 KMin = VCoord->MinLayerCellH(i);
         I4 KMax = std::min(VCoord->MaxLayerCellH(i), Dims[1] - 1);
         for (I4 k = KMin; k <= KMax; k++) {
            // Compute expected using same precision path as operator
            Real PseudoVal  = static_cast<Real>(InputHost(i, k));
            Real SpecVolVal = SpecVolHost(i, k);
            ScalarT Expected =
                static_cast<ScalarT>(RhoSw * SpecVolVal * PseudoVal);
            ScalarT Computed = ResultHost(i, k);

            if (std::abs(static_cast<Real>(Computed) -
                         static_cast<Real>(Expected)) /
                    static_cast<Real>(Expected) >
                Tolerance) {
               Passed = false;
               break;
            }
            NumVerified++;
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      // 3D: check all (i,j,k) points in active layers
      for (I4 i = 0; i < Dims[0]; i++) {
         for (I4 j = 0; j < Mesh->NCellsOwned; j++) {
            // Respect active layer bounds for this cell, clamped to array
            // bounds
            I4 KMin = VCoord->MinLayerCellH(j);
            I4 KMax = std::min(VCoord->MaxLayerCellH(j), Dims[2] - 1);
            for (I4 k = KMin; k <= KMax; k++) {
               // Compute expected using same precision path as operator
               Real PseudoVal  = static_cast<Real>(InputHost(i, j, k));
               Real SpecVolVal = SpecVolHost(j, k);
               ScalarT Expected =
                   static_cast<ScalarT>(RhoSw * SpecVolVal * PseudoVal);
               ScalarT Computed = ResultHost(i, j, k);

               if (std::abs(static_cast<Real>(Computed) -
                            static_cast<Real>(Expected)) /
                       static_cast<Real>(Expected) >
                   Tolerance) {
                  Passed = false;
                  break;
               }
               NumVerified++;
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   reportTest("PseudoToGeometricOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("PseudoToGeometricOp {} test failed after verifying {} points",
                TypeName, NumVerified);
   }
}

//------------------------------------------------------------------------------
// Test PseudoToGeometricOp with all supported array types
void testPseudoToGeometricOp(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   // Test 1D arrays
   testPseudoToGeometricOpType<Array1DR4>("Array1DR4", Env, Mesh, VCoord);
   testPseudoToGeometricOpType<Array1DR8>("Array1DR8", Env, Mesh, VCoord);

   // Test 2D arrays
   testPseudoToGeometricOpType<Array2DR4>("Array2DR4", Env, Mesh, VCoord);
   testPseudoToGeometricOpType<Array2DR8>("Array2DR8", Env, Mesh, VCoord);

   // Test 3D arrays
   testPseudoToGeometricOpType<Array3DR4>("Array3DR4", Env, Mesh, VCoord);
   testPseudoToGeometricOpType<Array3DR8>("Array3DR8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test ExtractRegionOp for a specific array type
template <typename ArrayType>
void testExtractRegionOpType(const std::string &TypeName, const MachEnv *Env,
                             const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   // Get dimensions and create test input field
   auto Dims     = Helper::getDims(Mesh, VCoord);
   auto DimNames = Helper::getDimNames();

   std::string InputFieldName = "TestExtractRegionInput_" + TypeName;
   std::string MaskFieldName  = "TestExtractRegionMask";

   // Create input field with known values
   Helper::createField(InputFieldName, Dims, [](auto... indices) {
      // Sum of indices + 1 to ensure non-zero values
      return static_cast<ScalarT>((indices + ...) + 1);
   });

   // Create mask field (1D horizontal, I4)
   // Mask = 1 for first half of cells, 0 for second half
   auto MaskField = Field::create(MaskFieldName, "Test regional mask", "1", "",
                                  0, 1, 1, {"NCells"});

   Array1DI4 MaskData(MaskFieldName + "_data", Mesh->NCellsSize);
   MaskField->attachData<Array1DI4>(MaskData, false);

   auto MaskDataHost = Kokkos::create_mirror_view(MaskData);
   I4 MidPoint       = Mesh->NCellsOwned / 2;
   for (I4 i = 0; i < Mesh->NCellsSize; ++i) {
      MaskDataHost(i) = (i < MidPoint) ? 1 : 0;
   }
   deepCopy(MaskData, MaskDataHost);

   // Create ExtractRegionOp
   auto Op = AnalysisOpFactory::createOp(
       "ExtractRegion", {InputFieldName},
       makeOpConfig(opParam("MaskName", MaskFieldName)));

   Op->initialize(Env, Mesh, VCoord, makeOpConfig());

   // Compute
   TimeInstant TestTime;
   Op->compute(TestTime);

   // Get result
   auto OutputNames = Op->getOutputFieldNames();
   OMEGA_ASSERT(OutputNames.size() == 1,
                "ExtractRegionOp should have 1 output");

   auto ResultField = Field::get(OutputNames[0]);
   auto ResultData  = ResultField->getDataArray<ArrayType>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Get input for verification
   auto InputField = Field::get(InputFieldName);
   auto InputData  = InputField->getDataArray<ArrayType>();
   auto InputHost  = createHostMirrorCopy(InputData);

   // Verify: where mask=1, result should equal input; where mask=0, result
   // should equal fill value
   bool Passed               = true;
   ScalarT Tol               = Helper::getTolerance();
   int NumVerified           = 0;
   constexpr ScalarT FillVal = FillValue<ScalarT>;

   if constexpr (Rank == 1) {
      // 1D: check all horizontal points
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         I4 MaskVal       = MaskDataHost(i);
         ScalarT Computed = ResultHost(i);
         ScalarT Expected;

         if (MaskVal == 0) {
            // Where mask is 0, expect fill value (not computed)
            Expected = FillVal;
         } else {
            // Where mask is 1, expect input value
            Expected = InputHost(i);
         }

         Real RelErr = std::abs(static_cast<Real>(Computed - Expected)) /
                       std::abs(static_cast<Real>(Expected));
         if (RelErr > Tol) {
            Passed = false;
            LOG_ERROR("  1D mismatch at cell {}: Mask={}, Expected {}, Got {}, "
                      "RelErr={}",
                      i, MaskVal, static_cast<Real>(Expected),
                      static_cast<Real>(Computed), RelErr);
            break;
         }
         NumVerified++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: check all (i,k) points in active layers
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         I4 KMin    = VCoord->MinLayerCellH(i);
         I4 KMax    = VCoord->MaxLayerCellH(i);
         I4 MaskVal = MaskDataHost(i);
         for (I4 k = KMin; k <= KMax; ++k) {
            ScalarT Computed = ResultHost(i, k);
            ScalarT Expected;

            if (MaskVal == 0) {
               // Where mask is 0, expect fill value (not computed)
               Expected = FillVal;
            } else {
               // Where mask is 1, expect input value
               Expected = InputHost(i, k);
            }

            Real RelErr = std::abs(static_cast<Real>(Computed - Expected)) /
                          std::abs(static_cast<Real>(Expected));
            if (RelErr > Tol) {
               Passed = false;
               LOG_ERROR("  2D mismatch at cell {}, layer {}: Mask={}, "
                         "Expected {}, Got {}, RelErr={}",
                         i, k, MaskVal, static_cast<Real>(Expected),
                         static_cast<Real>(Computed), RelErr);
               break;
            }
            NumVerified++;
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      // 3D: check all (i,j,k) points in active layers
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            I4 KMin    = VCoord->MinLayerCellH(j);
            I4 KMax    = VCoord->MaxLayerCellH(j);
            I4 MaskVal = MaskDataHost(j);
            for (I4 k = KMin; k <= KMax; ++k) {
               ScalarT Computed = ResultHost(i, j, k);
               ScalarT Expected;

               if (MaskVal == 0) {
                  // Where mask is 0, expect fill value (not computed)
                  Expected = FillVal;
               } else {
                  // Where mask is 1, expect input value
                  Expected = InputHost(i, j, k);
               }

               Real RelErr = std::abs(static_cast<Real>(Computed - Expected)) /
                             std::abs(static_cast<Real>(Expected));
               if (RelErr > Tol) {
                  Passed = false;
                  LOG_ERROR("  3D mismatch at i={}, cell {}, layer {}: "
                            "Mask={}, Expected {}, Got {}, RelErr={}",
                            i, j, k, MaskVal, static_cast<Real>(Expected),
                            static_cast<Real>(Computed), RelErr);
                  break;
               }
               NumVerified++;
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   // Verify output field has mask attached
   if (Passed && !ResultField->hasRegionalMask()) {
      Passed = false;
      LOG_ERROR(
          "  ExtractRegionOp output field missing regional mask attachment");
   }

   reportTest("ExtractRegionOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("ExtractRegionOp {} test failed after verifying {} points",
                TypeName, NumVerified);
   }

   Field::destroy(MaskFieldName);
}

//------------------------------------------------------------------------------
// Test ExtractRegionOp with all supported array types
void testExtractRegionOp(const MachEnv *Env, const HorzMesh *Mesh,
                         const VertCoord *VCoord) {

   // Test 1D arrays (all scalar types)
   testExtractRegionOpType<Array1DI4>("Array1DI4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array1DI8>("Array1DI8", Env, Mesh, VCoord);
   testExtractRegionOpType<Array1DR4>("Array1DR4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array1DR8>("Array1DR8", Env, Mesh, VCoord);

   // Test 2D arrays (all scalar types)
   testExtractRegionOpType<Array2DI4>("Array2DI4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array2DI8>("Array2DI8", Env, Mesh, VCoord);
   testExtractRegionOpType<Array2DR4>("Array2DR4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array2DR8>("Array2DR8", Env, Mesh, VCoord);

   // Test 3D arrays (all scalar types)
   testExtractRegionOpType<Array3DI4>("Array3DI4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array3DI8>("Array3DI8", Env, Mesh, VCoord);
   testExtractRegionOpType<Array3DR4>("Array3DR4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array3DR8>("Array3DR8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test TransectAccumulatorOp with edge-based transport
// Tests the accumulation of edge transport across a synthetic transect.
// Uses non-contiguous edge pattern (every 5th edge) to thoroughly test
// MPI reduction; real transects would be spatially contiguous.
void testTransectAccumulatorOp(const MachEnv *Env, const HorzMesh *Mesh,
                               const VertCoord *VCoord) {

   I4 NEdges      = Mesh->NEdgesSize;
   I4 NVertLayers = VCoord->NVertLayers;

   // Create transport field (NEdges × NVertLayers)
   // Value at each point = layer index + 1 (so layer 0→1.0, layer 1→2.0, etc.)
   std::string TransportFieldName = "TestTransectTransport";
   auto TransportField =
       Field::create(TransportFieldName, "Test transport for transect", "m3/s",
                     "", -1e30, 1e30, 2, {"NEdges", "NVertLayers"});

   Array2DR8 TransportData(TransportFieldName + "_data", NEdges, NVertLayers);
   TransportField->attachData<Array2DR8>(TransportData, false);

   auto TransportDataHost = Kokkos::create_mirror_view(TransportData);
   for (I4 i = 0; i < NEdges; ++i) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         TransportDataHost(i, k) =
             static_cast<Real>(k + 1); // Layer-dependent value
      }
   }
   deepCopy(TransportData, TransportDataHost);

   // Create transect mask field (which edges are on transect)
   // Every 5th edge is on transect (mask=1.0), others are off (mask=0.0)
   std::string MaskFieldName = "TestTransectMask";
   auto MaskField = Field::create(MaskFieldName, "Test transect mask", "1", "",
                                  0, 1, 1, {"NEdges"});

   Array1DI4 MaskData(MaskFieldName + "_data", NEdges);
   MaskField->attachData<Array1DI4>(MaskData, false);

   auto MaskDataHost = Kokkos::create_mirror_view(MaskData);
   for (I4 i = 0; i < NEdges; ++i) {
      MaskDataHost(i) = ((i % 5) == 0) ? 1 : 0;
   }
   deepCopy(MaskData, MaskDataHost);

   // Create transect sign field (direction of transport)
   // Alternating +1.0 and -1.0 for edges on transect
   std::string SignFieldName = "TestTransectSign";
   auto SignField = Field::create(SignFieldName, "Test transect sign", "1", "",
                                  -1, 1, 1, {"NEdges"});

   Array1DI4 SignData(SignFieldName + "_data", NEdges);
   SignField->attachData<Array1DI4>(SignData, false);

   auto SignDataHost = Kokkos::create_mirror_view(SignData);
   for (I4 i = 0; i < NEdges; ++i) {
      if ((i % 5) == 0) {
         // Edges on transect: alternate sign based on edge index
         SignDataHost(i) = ((i / 5) % 2 == 0) ? 1 : -1;
      } else {
         SignDataHost(i) = 0; // Doesn't matter, mask is 0
      }
   }
   deepCopy(SignData, SignDataHost);

   // Calculate expected result per layer
   // For each layer k: sum over owned edges of transport(i,k) × mask(i) ×
   // sign(i)
   std::vector<Real> LocalSumPerLayer(NVertLayers, 0.0);

   for (I4 i = 0; i < Mesh->NEdgesOwned; ++i) {
      if (MaskDataHost(i) != 0) { // Edge is on transect
         I4 KMin = VCoord->MinLayerEdgeBotH(i);
         I4 KMax = VCoord->MaxLayerEdgeTopH(i);
         for (I4 k = KMin; k <= KMax; ++k) {
            LocalSumPerLayer[k] +=
                TransportDataHost(i, k) * static_cast<Real>(SignDataHost(i));
         }
      }
   }

   // Global sum across MPI ranks
   std::vector<Real> ExpectedPerLayer(NVertLayers);
   for (I4 k = 0; k < NVertLayers; ++k) {
      ExpectedPerLayer[k] = globalSum(LocalSumPerLayer[k], Env->getComm());
   }

   // Create TransectAccumulator operator
   auto TransectOp = AnalysisOpFactory::createOp(
       "TransectAccumulator",
       {TransportFieldName, MaskFieldName, SignFieldName},
       makeOpConfig(opParam("TransectName", "TestTransect")));
   TransectOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   TransectOp->compute(TestTime);

   // Get result (should be 1D: NVertLayers)
   // Use operator's OutputNames to get the correct field name
   auto OutputNames = TransectOp->getOutputFieldNames();
   if (OutputNames.empty()) {
      LOG_ERROR("TransectAccumulatorOp has no output names!");
      reportTest("TransectAccumulatorOp", false);
      return;
   }
   auto ResultField = Field::get(OutputNames[0]);
   auto ResultData  = ResultField->getDataArray<Array1DR8>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Verify results for all layers
   bool Passed = true;
   for (I4 k = 0; k < NVertLayers; ++k) {
      Real Computed = ResultHost(k);
      Real Expected = ExpectedPerLayer[k];

      if (std::abs(Computed - Expected) > 1.0e-8) {
         Passed = false;
         LOG_ERROR("  Layer {}: Expected {}, Got {}", k, Expected, Computed);
         break;
      }
   }

   reportTest("TransectAccumulatorOp", Passed);
}

//===----------------------------------------------------------------------===//
// Initialization and finalization functions
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Initialize needed modules
void initAnalysisTest() {

   I4 Err;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();

   // Get the model clock
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize streams
   IOStream::init(ModelClock);

   // Initialize fields
   Field::init(ModelClock);

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0)
      ABORT_ERROR("AnalysisOperatorTest: error initializing default halo");

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   VertCoord::init();

   // Initialize tracers
   Tracers::init();

   // Initialize auxiliary state
   AuxiliaryState::init();

   // Initialize equation of state
   Eos::init();

   // Initialize pressure gradient
   PressureGrad::init();

   // Initialize forcing
   Forcing::init();

   // Initialize vertical mixing
   VertMix::init();

   // Initialize tendencies
   Tendencies::init();

   // Initialize vertical advection
   VertAdv::init();

   // Second step of time stepper initialization
   TimeStepper::init2();

   // Initialize ocean state
   Err = OceanState::init();
   if (Err != 0)
      ABORT_ERROR("AnalysisOperatorTest: error initializing default state");

   // Register all analysis operators
   Analysis::init();
}

//------------------------------------------------------------------------------
// Clean-up modules
void finalizeAnalysisTest() {

   Analysis::finalize();
   IOStream::finalize();
   VertMix::destroyInstance();
   Forcing::clear();
   OceanState::clear();
   Tracers::clear();
   AuxiliaryState::clear();
   PressureGrad::clear();
   Tendencies::clear();
   VertAdv::clear();
   VertCoord::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

//===----------------------------------------------------------------------===//
// Main test driver
//===----------------------------------------------------------------------===//

int main(int argc, char *argv[]) {

   int Err = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      initAnalysisTest();

      auto DefEnv     = MachEnv::getDefault();
      auto DefStepper = TimeStepper::getDefault();
      auto Mesh       = HorzMesh::getDefault();
      auto VCoord     = VertCoord::getDefault();
      auto ModelClock = DefStepper->getClock();

      testSpatialMaxOp(DefEnv, Mesh, VCoord);

      testSpatialMinOp(DefEnv, Mesh, VCoord);

      testSpatialMeanOp(DefEnv, Mesh, VCoord);

      testSpatialStdDevOp(DefEnv, Mesh, VCoord);

      testTimeMeanOp(DefEnv, Mesh, VCoord, ModelClock);

      testScalarMultiplyOp(DefEnv, Mesh, VCoord);

      testBinaryMultiplyOp(DefEnv, Mesh, VCoord);

      testPrefixSumOp(DefEnv, Mesh, VCoord);

      testCoordinateBinningOp(DefEnv, Mesh, VCoord);

      testBinnedAccumulatorOp(DefEnv, Mesh, VCoord);

      testPseudoToGeometricOp(DefEnv, Mesh, VCoord);

      testExtractRegionOp(DefEnv, Mesh, VCoord);

      testTransectAccumulatorOp(DefEnv, Mesh, VCoord);

      if (NumFailed > 0) {
         Err = 1;
         LOG_ERROR("AnalysisOperatorTest failure");
         LOG_ERROR("  Total tests: {}", NumTests);
         LOG_ERROR("  Passed: {}", NumPassed);
         LOG_ERROR("  Failed: {}", NumFailed);
      }

      finalizeAnalysisTest();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return Err;
}
