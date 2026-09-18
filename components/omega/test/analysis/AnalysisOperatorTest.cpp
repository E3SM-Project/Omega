//===-- Test 5.1: Multi-type operator correctness tests --------*- C++ -*-===//
//
// Comprehensive unit tests for Analysis operators across all supported
// array types (ranks 1D/2D/3D and scalar types I4/I8/R4/R8)
//
//===-----------------------------------------------------------------------===//

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
#include "OceanState.h"
#include "TimeStepper.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "VertMix.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
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

   // Relative tolerance for a weighted statistic, whose summation order
   // differs from the host-side check
   static Real getRelTolerance() {
      if constexpr (std::is_same_v<ScalarT, float>) {
         return 1.0e-5;
      } else {
         return 1.0e-12;
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
      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", -1.0e30, 1.0e30, 1, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0]);
      TestField->attachData<ArrayType>(TestData);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         TestDataHost(i) = ValueFunc(i);
      }
      Kokkos::deep_copy(TestData, TestDataHost);
   }

   // Create test field for 2D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 2, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();
      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", -1.0e30, 1.0e30, 2, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
      TestField->attachData<ArrayType>(TestData);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            TestDataHost(i, j) = ValueFunc(i, j);
         }
      }
      Kokkos::deep_copy(TestData, TestDataHost);
   }

   // Create test field for 3D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 3, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();
      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", -1.0e30, 1.0e30, 3, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1], Dims[2]);
      TestField->attachData<ArrayType>(TestData);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            for (I4 k = 0; k < Dims[2]; ++k) {
               TestDataHost(i, j, k) = ValueFunc(i, j, k);
            }
         }
      }
      Kokkos::deep_copy(TestData, TestDataHost);
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

//------------------------------------------------------------------------------
// Host-side weights for checking the spatial statistics, computed from the
// host mesh, mask and pseudo-thickness arrays independently of SpatialWeights:
// area for horizontal fields, mass (area times pseudo-thickness) for layered
// fields with the thickness interpolated to edges and vertices as the
// auxiliary variables do, and half the mass of each adjacent layer for
// interface fields. Inactive entries have zero weight.
struct HostWeights {
   const HorzMesh *Mesh;
   const VertCoord *VCoord;
   HostArray2DReal ThickH;

   HostWeights(const HorzMesh *InMesh, const VertCoord *InVCoord)
       : Mesh(InMesh), VCoord(InVCoord) {
      auto State      = OceanState::getDefault();
      auto ThickField = Field::get(State->PseudoThicknessFldName);
      ThickH = createHostMirrorCopy(ThickField->getDataArray<Array2DReal>());
   }

   bool cellActive(I4 ICell, I4 K) const {
      return VCoord->CellMaskH(ICell, K) > 0;
   }

   // Mass of one layer of one entity
   Real cellLayer(I4 ICell, I4 K) const {
      return cellActive(ICell, K) ? Mesh->AreaCellH(ICell) * ThickH(ICell, K)
                                  : 0;
   }
   Real edgeLayer(I4 IEdge, I4 K) const {
      if (VCoord->EdgeMaskH(IEdge, K) <= 0)
         return 0;
      Real Thick = 0;
      for (I4 J = 0; J < 2; ++J) {
         I4 JCell = Mesh->CellsOnEdgeH(IEdge, J);
         if (cellActive(JCell, K))
            Thick += 0.5 * ThickH(JCell, K);
      }
      return Mesh->DcEdgeH(IEdge) * Mesh->DvEdgeH(IEdge) * Thick;
   }
   Real vertexLayer(I4 IVertex, I4 K) const {
      if (VCoord->VertexMaskH(IVertex, K) <= 0)
         return 0;
      Real Mass = 0;
      for (I4 J = 0; J < Mesh->VertexDegree; ++J) {
         I4 JCell = Mesh->CellsOnVertexH(IVertex, J);
         if (cellActive(JCell, K))
            Mass += Mesh->KiteAreasOnVertexH(IVertex, J) * ThickH(JCell, K);
      }
      return Mass;
   }

   // Weight of one interface of one cell: half of each adjacent layer
   Real cellInterface(I4 ICell, I4 K) const {
      Real Above = (K > 0) ? cellLayer(ICell, K - 1) : 0;
      Real Below = (K < VCoord->NVertLayers) ? cellLayer(ICell, K) : 0;
      return 0.5 * (Above + Below);
   }

   // Area weights of horizontal fields
   Real cellArea(I4 ICell) const {
      return cellActive(ICell, 0) ? Mesh->AreaCellH(ICell) : 0;
   }
   Real edgeArea(I4 IEdge) const {
      return (VCoord->EdgeMaskH(IEdge, 0) > 0)
                 ? Mesh->DcEdgeH(IEdge) * Mesh->DvEdgeH(IEdge)
                 : 0;
   }
};

//------------------------------------------------------------------------------
// Computes the expected weighted mean and standard deviation of a field from
// a visitor that calls back with the weight and value of every local entry.
// Two passes, so the deviations are taken from the global mean.
using WeightedVisitor = std::function<void(std::function<void(Real, Real)>)>;

void expectedStats(const WeightedVisitor &Visit, MPI_Comm Comm, Real &Mean,
                   Real &StdDev) {
   Real SumW = 0, SumWX = 0;
   Visit([&](Real W, Real X) {
      if (W > 0) {
         SumW += W;
         SumWX += W * X;
      }
   });
   Real W = globalSum(SumW, Comm);
   Mean   = globalSum(SumWX, Comm) / W;

   Real SumWDev = 0;
   Visit([&](Real Wt, Real X) {
      if (Wt > 0)
         SumWDev += Wt * (X - Mean) * (X - Mean);
   });
   StdDev = std::sqrt(globalSum(SumWDev, Comm) / W);
}

//------------------------------------------------------------------------------
// Checks a computed statistic against its expected value to a relative
// tolerance
bool checkClose(const std::string &TestName, Real Computed, Real Expected,
                Real RelTol) {
   Real Scale  = std::max(std::abs(Expected), static_cast<Real>(1));
   bool Passed = std::abs(Computed - Expected) <= RelTol * Scale;
   reportTest(TestName, Passed);
   if (!Passed)
      LOG_ERROR("  Expected: {}, Got: {}", Expected, Computed);
   return Passed;
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
   Config EmptyConfig;
   auto MaxOp =
       AnalysisOpFactory::createOp("SpatialMax", {FieldName}, EmptyConfig);
   MaxOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MaxOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMax");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

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
   Config EmptyConfig;
   auto MinOp =
       AnalysisOpFactory::createOp("SpatialMin", {FieldName}, EmptyConfig);
   MinOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MinOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMin");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

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

   // Expected weighted mean from the host-side weights: area for the 1D
   // field, mass for the layered ones, over owned cells and active layers
   HostWeights Weights(Mesh, VCoord);
   auto Value = [Value1, Value2](I4 Parity) -> Real {
      return static_cast<Real>((Parity % 2) == 0 ? Value1 : Value2);
   };
   WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
      if constexpr (Rank == 1) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            Add(Weights.cellArea(i), Value(i));
      } else if constexpr (Rank == 2) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            for (I4 j = 0; j < VCoord->NVertLayers; ++j)
               Add(Weights.cellLayer(i, j), Value(i + j));
      } else if constexpr (Rank == 3) {
         I4 NTracers = Tracers::getNumTracers();
         for (I4 t = 0; t < NTracers; ++t)
            for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
               for (I4 j = 0; j < VCoord->NVertLayers; ++j)
                  Add(Weights.cellLayer(i, j), Value(t + i + j));
      }
   };
   Real ExpectedMean, ExpectedStdDev;
   expectedStats(Visit, Env->getComm(), ExpectedMean, ExpectedStdDev);

   // Create and compute operator
   Config EmptyConfig;
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialMean");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

   Real ComputedMean = ResultHost(0);

   // Verify
   checkClose("SpatialMeanOp: " + TypeName, ComputedMean, ExpectedMean,
              Helper::getRelTolerance());
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

   // Expected weighted standard deviation from the host-side weights
   HostWeights Weights(Mesh, VCoord);
   auto Value = [Value1, Value2](I4 Parity) -> Real {
      return static_cast<Real>((Parity % 2) == 0 ? Value1 : Value2);
   };
   WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
      if constexpr (Rank == 1) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            Add(Weights.cellArea(i), Value(i));
      } else if constexpr (Rank == 2) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            for (I4 j = 0; j < VCoord->NVertLayers; ++j)
               Add(Weights.cellLayer(i, j), Value(i + j));
      } else if constexpr (Rank == 3) {
         I4 NTracers = Tracers::getNumTracers();
         for (I4 t = 0; t < NTracers; ++t)
            for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
               for (I4 j = 0; j < VCoord->NVertLayers; ++j)
                  Add(Weights.cellLayer(i, j), Value(t + i + j));
      }
   };
   Real ExpectedMean, ExpectedStdDev;
   expectedStats(Visit, Env->getComm(), ExpectedMean, ExpectedStdDev);

   // SpatialStdDevOp requires a pre-existing _SpatialMean field for the input.
   // Create and compute a SpatialMeanOp first so that field is registered.
   Config EmptyConfig;
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Now create and compute the StdDev operator
   auto StdDevOp =
       AnalysisOpFactory::createOp("SpatialStdDev", {FieldName}, EmptyConfig);
   StdDevOp->initialize(Env, Mesh, VCoord, EmptyConfig);
   StdDevOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialStdDev");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

   Real ComputedStdDev = ResultHost(0);

   // Verify
   checkClose("SpatialStdDevOp: " + TypeName, ComputedStdDev, ExpectedStdDev,
              Helper::getRelTolerance());
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
      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      Kokkos::deep_copy(TestDataHost, TestData);

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

      Kokkos::deep_copy(TestData, TestDataHost);

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
      auto ResultHost = Kokkos::create_mirror_view(ResultData);
      Kokkos::deep_copy(ResultHost, ResultData);

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
      auto ResultHost = Kokkos::create_mirror_view(ResultData);
      Kokkos::deep_copy(ResultHost, ResultData);

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
      auto ResultHost = Kokkos::create_mirror_view(ResultData);
      Kokkos::deep_copy(ResultHost, ResultData);

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

//------------------------------------------------------------------------------
// Creates a 2D Real field on the given entities with the given vertical
// dimension, filled with a pattern of values on active owned entries and
// NaN everywhere else (inactive layers and halo entities), so that a
// statistic that reads an inactive or halo entry comes out NaN. Returns the
// host copy of the values.
HostArray2DReal createEntityField(const std::string &FieldName,
                                  const std::string &HorizDim,
                                  const std::string &VertDim, I4 NEntities,
                                  I4 NOwned, I4 NVert,
                                  std::function<bool(I4, I4)> IsActive,
                                  std::function<Real(I4, I4)> ValueFunc) {

   std::vector<std::string> DimNames = {HorizDim, VertDim};
   auto TestField = Field::create(FieldName, "Weighted statistics test field",
                                  "m", "", -1.0e30, 1.0e30, 2, DimNames);
   Array2DReal Data(FieldName + "_data", NEntities, NVert);
   TestField->attachData<Array2DReal>(Data);

   auto DataHost  = Kokkos::create_mirror_view(Data);
   const Real NaN = std::numeric_limits<Real>::quiet_NaN();
   for (I4 I = 0; I < NEntities; ++I) {
      for (I4 K = 0; K < NVert; ++K) {
         bool Use       = (I < NOwned) && IsActive(I, K);
         DataHost(I, K) = Use ? ValueFunc(I, K) : NaN;
      }
   }
   Kokkos::deep_copy(Data, DataHost);
   return DataHost;
}

//------------------------------------------------------------------------------
// Computes the mean and standard deviation of a field with the operators and
// checks them against the expected values
void checkMeanAndStdDev(const std::string &TestName,
                        const std::string &FieldName, const MachEnv *Env,
                        const HorzMesh *Mesh, const VertCoord *VCoord,
                        Real ExpectedMean, Real ExpectedStdDev, Real RelTol) {

   Config EmptyConfig;
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);
   auto StdDevOp =
       AnalysisOpFactory::createOp("SpatialStdDev", {FieldName}, EmptyConfig);
   StdDevOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MeanOp->compute(TestTime);
   StdDevOp->compute(TestTime);

   auto MeanHost = createHostMirrorCopy(
       Field::get(FieldName + "_SpatialMean")->getDataArray<Array1DReal>());
   auto StdDevHost = createHostMirrorCopy(
       Field::get(FieldName + "_SpatialStdDev")->getDataArray<Array1DReal>());

   checkClose(TestName + " mean", MeanHost(0), ExpectedMean, RelTol);
   checkClose(TestName + " std dev", StdDevHost(0), ExpectedStdDev, RelTol);
}

//------------------------------------------------------------------------------
// Tests the weighted mean and standard deviation on edges, vertices and cell
// interfaces, and on a horizontal edge field, against host-side weights.
// Inactive and halo entries hold NaN, so a statistic that touched them would
// fail. A constant field must give its value and a zero standard deviation
// whatever the weights.
void testWeightedStats(const MachEnv *Env, const HorzMesh *Mesh,
                       const VertCoord *VCoord) {

   HostWeights Weights(Mesh, VCoord);
   MPI_Comm Comm      = Env->getComm();
   const Real RelTol  = 1.0e-12;
   const I4 NVertLyrs = VCoord->NVertLayers;
   auto Pattern       = [](I4 I, I4 K) -> Real {
      return static_cast<Real>(10 + (7 * I + 3 * K) % 5);
   };
   auto Constant = [](I4 I, I4 K) -> Real { return static_cast<Real>(7.5); };
   Real ExpectedMean, ExpectedStdDev;

   // Layered field on edges, mass-weighted with the edge thickness
   {
      auto Active = [&](I4 E, I4 K) { return VCoord->EdgeMaskH(E, K) > 0; };
      createEntityField("TestWeightedEdge", "NEdges", "NVertLayers",
                        Mesh->NEdgesSize, Mesh->NEdgesOwned, NVertLyrs, Active,
                        Pattern);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 E = 0; E < Mesh->NEdgesOwned; ++E)
            for (I4 K = 0; K < NVertLyrs; ++K)
               Add(Weights.edgeLayer(E, K), Pattern(E, K));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: edge layers", "TestWeightedEdge", Env,
                         Mesh, VCoord, ExpectedMean, ExpectedStdDev, RelTol);
   }

   // Layered field on vertices, mass-weighted with the kite areas
   {
      auto Active = [&](I4 V, I4 K) { return VCoord->VertexMaskH(V, K) > 0; };
      createEntityField("TestWeightedVertex", "NVertices", "NVertLayers",
                        Mesh->NVerticesSize, Mesh->NVerticesOwned, NVertLyrs,
                        Active, Pattern);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 V = 0; V < Mesh->NVerticesOwned; ++V)
            for (I4 K = 0; K < NVertLyrs; ++K)
               Add(Weights.vertexLayer(V, K), Pattern(V, K));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: vertex layers", "TestWeightedVertex",
                         Env, Mesh, VCoord, ExpectedMean, ExpectedStdDev,
                         RelTol);
   }

   // Interface field on cells, half the mass of each adjacent layer. The
   // bottom interface of every active column takes part.
   {
      auto Active = [&](I4 C, I4 K) {
         bool Above = (K > 0) && Weights.cellActive(C, K - 1);
         bool Below = (K < NVertLyrs) && Weights.cellActive(C, K);
         return Above || Below;
      };
      createEntityField("TestWeightedInterface", "NCells", "NVertLayersP1",
                        Mesh->NCellsSize, Mesh->NCellsOwned, NVertLyrs + 1,
                        Active, Pattern);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 C = 0; C < Mesh->NCellsOwned; ++C)
            for (I4 K = 0; K <= NVertLyrs; ++K)
               Add(Weights.cellInterface(C, K), Pattern(C, K));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: cell interfaces",
                         "TestWeightedInterface", Env, Mesh, VCoord,
                         ExpectedMean, ExpectedStdDev, RelTol);
   }

   // Horizontal field on edges, area-weighted
   {
      std::vector<std::string> DimNames = {"NEdges"};
      auto TestField =
          Field::create("TestWeightedEdge1D", "Weighted statistics test field",
                        "m", "", -1.0e30, 1.0e30, 1, DimNames);
      Array1DReal Data("TestWeightedEdge1D_data", Mesh->NEdgesSize);
      TestField->attachData<Array1DReal>(Data);
      auto DataHost  = Kokkos::create_mirror_view(Data);
      const Real NaN = std::numeric_limits<Real>::quiet_NaN();
      for (I4 E = 0; E < Mesh->NEdgesSize; ++E) {
         bool Use    = (E < Mesh->NEdgesOwned) && (Weights.edgeArea(E) > 0);
         DataHost(E) = Use ? Pattern(E, 0) : NaN;
      }
      Kokkos::deep_copy(Data, DataHost);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 E = 0; E < Mesh->NEdgesOwned; ++E)
            Add(Weights.edgeArea(E), Pattern(E, 0));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: horizontal edges",
                         "TestWeightedEdge1D", Env, Mesh, VCoord, ExpectedMean,
                         ExpectedStdDev, RelTol);
   }

   // A constant layered field on cells: the mean is the constant and the
   // standard deviation zero, whatever the weights
   {
      auto Active = [&](I4 C, I4 K) { return Weights.cellActive(C, K); };
      createEntityField("TestWeightedConstant", "NCells", "NVertLayers",
                        Mesh->NCellsSize, Mesh->NCellsOwned, NVertLyrs, Active,
                        Constant);
      checkMeanAndStdDev("Weighted stats: constant field",
                         "TestWeightedConstant", Env, Mesh, VCoord, 7.5, 0.0,
                         RelTol);
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

   // Initialize VertAdv
   VertAdv::init();

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

   // Read the initial state so the pseudo-thickness the mass weights use is
   // the one from the mesh file rather than fill values
   if (!IOStream::validateAll())
      ABORT_ERROR("AnalysisOperatorTest: stream validation failed");
   Metadata ReqMeta;
   Error ReadErr = IOStream::read("InitialState", ModelClock, ReqMeta);
   CHECK_ERROR_ABORT(ReadErr,
                     "AnalysisOperatorTest: failed to read initial state");

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

      testWeightedStats(DefEnv, Mesh, VCoord);

      testTimeMeanOp(DefEnv, Mesh, VCoord, ModelClock);

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
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return Err;
}
