#include "HorzOperators.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <cmath>

using namespace OMEGA;

// check if two real numbers are equal with a given relative tolerance
bool isApprox(Real X, Real Y, Real RTol) {
   return std::abs(X - Y) <= RTol * std::max(std::abs(X), std::abs(Y));
}

// temporary replacement for YAKL intrinsics
Real maxVal(const Array1DReal &Arr) {
   Real MaxVal;

   parallelReduce(
       {Arr.extent_int(0)},
       KOKKOS_LAMBDA(int I, Real &Accum) {
          Accum = Kokkos::max(Arr(I), Accum);
       },
       Kokkos::Max<Real>(MaxVal));

   return MaxVal;
}

Real sum(const Array1DReal &Arr) {
   Real Sum;

   parallelReduce(
       {Arr.extent_int(0)},
       KOKKOS_LAMBDA(int I, Real &Accum) { Accum += Arr(I); }, Sum);

   return Sum;
}

#ifdef HORZOPERATORS_TEST_PLANE
// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetup {
   Real Pi = M_PI;

   // lengths of periodic planar mesh
   // TODO: get this from the horizontal mesh once it supports periodic planar
   // meshes
   Real Lx = 1;
   Real Ly = std::sqrt(3) / 2;

   Real ExpectedDivErrorLInf   = 0.00124886886594427027;
   Real ExpectedDivErrorL2     = 0.00124886886590974385;
   Real ExpectedGradErrorLInf  = 0.00125026071878537952;
   Real ExpectedGradErrorL2    = 0.00134354611117262204;
   Real ExpectedCurlErrorLInf  = 0.161365663569699946;
   Real ExpectedCurlErrorL2    = 0.161348016897141039;
   Real ExpectedReconErrorLInf = 0.00450897496974901352;
   Real ExpectedReconErrorL2   = 0.00417367308684470691;

   KOKKOS_INLINE_FUNCTION Real exactScalar(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_INLINE_FUNCTION Real exactGradScalarX(Real X, Real Y) const {
      return 2 * Pi / Lx * std::cos(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_INLINE_FUNCTION Real exactGradScalarY(Real X, Real Y) const {
      return 2 * Pi / Ly * std::sin(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_INLINE_FUNCTION Real exactVecX(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_INLINE_FUNCTION Real exactVecY(Real X, Real Y) const {
      return std::cos(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_INLINE_FUNCTION Real exactDivVec(Real X, Real Y) const {
      return 2 * Pi * (1. / Lx + 1. / Ly) * std::cos(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_INLINE_FUNCTION Real exactCurlVec(Real X, Real Y) const {
      return 2 * Pi * (-1. / Lx + 1. / Ly) * std::sin(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }
};
#endif

#ifdef HORZOPERATORS_TEST_SPHERE_1
// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetup {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = 6371220;

   Real ExpectedDivErrorLInf   = 0.0158660625383131929;
   Real ExpectedDivErrorL2     = 0.00401711826857833413;
   Real ExpectedGradErrorLInf  = 0.00182941720362535579;
   Real ExpectedGradErrorL2    = 0.00152823722290378484;
   Real ExpectedCurlErrorLInf  = 0.0361051145615213023;
   Real ExpectedCurlErrorL2    = 0.0271253226603560341;
   Real ExpectedReconErrorLInf = 0.0207855998352246864;
   Real ExpectedReconErrorL2   = 0.00687944381487612909;

   KOKKOS_INLINE_FUNCTION Real exactScalar(Real Lon, Real Lat) const {
      return Radius * std::cos(Lon) * std::pow(std::cos(Lat), 4);
   }

   KOKKOS_INLINE_FUNCTION Real exactGradScalarX(Real Lon, Real Lat) const {
      return -std::sin(Lon) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_INLINE_FUNCTION Real exactGradScalarY(Real Lon, Real Lat) const {
      return -4 * std::cos(Lon) * std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_INLINE_FUNCTION Real exactVecX(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_INLINE_FUNCTION Real exactVecY(Real Lon, Real Lat) const {
      return -4 * Radius * std::sin(Lon) * std::cos(Lon) *
             std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_INLINE_FUNCTION Real exactDivVec(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
             (20 * std::pow(std::sin(Lat), 2) - 6);
   }

   KOKKOS_INLINE_FUNCTION Real exactCurlVec(Real Lon, Real Lat) const {
      return -4 * std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 2) *
             std::sin(Lat);
   }
};
#endif

#ifdef HORZOPERATORS_TEST_SPHERE_2
// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetup {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = 6371220;

   Real ExpectedDivErrorLInf   = 1.45295480324903797e-09;
   Real ExpectedDivErrorL2     = 0.00429375897154735467;
   Real ExpectedGradErrorLInf  = 0.00202094598872142039;
   Real ExpectedGradErrorL2    = 0.00117807041607626765;
   Real ExpectedCurlErrorLInf  = 0.0328084017109396969;
   Real ExpectedCurlErrorL2    = 0.0105749467007983152;
   Real ExpectedReconErrorLInf = 0.0253216806569417016;
   Real ExpectedReconErrorL2   = 0.00425161856853827763;

   KOKKOS_INLINE_FUNCTION Real exactScalar(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lat), 2);
   }

   KOKKOS_INLINE_FUNCTION Real exactGradScalarX(Real Lon, Real Lat) const {
      return 0;
   }

   KOKKOS_INLINE_FUNCTION Real exactGradScalarY(Real Lon, Real Lat) const {
      return -2 * std::sin(Lat) * std::cos(Lat);
   }

   KOKKOS_INLINE_FUNCTION Real exactVecX(Real Lon, Real Lat) const {
      return std::cos(Lat);
   }

   KOKKOS_INLINE_FUNCTION Real exactVecY(Real Lon, Real Lat) const { return 0; }

   KOKKOS_INLINE_FUNCTION Real exactDivVec(Real Lon, Real Lat) const {
      return 0;
   }

   KOKKOS_INLINE_FUNCTION Real exactCurlVec(Real Lon, Real Lat) const {
      return 2 * std::sin(Lat) / Radius;
   }
};
#endif

int testDivergence(Real RTol) {
   int Err;
   TestSetup Setup;

   const auto &mesh = HorzMesh::getDefault();
#ifdef HORZOPERATORS_TEST_PLANE
   auto XEdge = createDeviceMirrorCopy(mesh->XEdgeH);
   auto YEdge = createDeviceMirrorCopy(mesh->YEdgeH);
#else
   auto XEdge = createDeviceMirrorCopy(mesh->LonEdgeH);
   auto YEdge = createDeviceMirrorCopy(mesh->LatEdgeH);
#endif
   auto &AngleEdge = mesh->AngleEdge;

   // Prepare operator input
   Array1DReal VecEdge("VecEdge", mesh->NEdgesSize);
   parallelFor(
       {mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          const Real X = XEdge(IEdge);
          const Real Y = YEdge(IEdge);

          const Real VecX = Setup.exactVecX(X, Y);
          const Real VecY = Setup.exactVecY(X, Y);

          const Real EdgeNormalX = std::cos(AngleEdge(IEdge));
          const Real EdgeNormalY = std::sin(AngleEdge(IEdge));

          VecEdge(IEdge) = EdgeNormalX * VecX + EdgeNormalY * VecY;
       });

   // Perform halo exchange
   Halo MyHalo(MachEnv::getDefaultEnv(), Decomp::getDefault());
   auto VecEdgeH = createHostMirrorCopy(VecEdge);
   MyHalo.exchangeFullArrayHalo(VecEdgeH, OnEdge);
   deepCopy(VecEdge, VecEdgeH);

#ifdef HORZOPERATORS_TEST_PLANE
   auto XCell = createDeviceMirrorCopy(mesh->XCellH);
   auto YCell = createDeviceMirrorCopy(mesh->YCellH);
#else
   auto XCell = createDeviceMirrorCopy(mesh->LonCellH);
   auto YCell = createDeviceMirrorCopy(mesh->LatCellH);
#endif
   auto &AreaCell = mesh->AreaCell;

   // Compute element-wise errors
   Array1DReal LInfCell("LInfCell", mesh->NCellsOwned);
   Array1DReal L2Cell("L2Cell", mesh->NCellsOwned);

   Array1DReal LInfScaleCell("LInfScaleCell", mesh->NCellsOwned);
   Array1DReal L2ScaleCell("L2ScaleCell", mesh->NCellsOwned);
   DivergenceOnCell DivergenceCell(mesh);
   parallelFor(
       {mesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          // Numerical result
          const Real DivCellNum = DivergenceCell(ICell, VecEdge);

          // Exact result
          const Real X            = XCell(ICell);
          const Real Y            = YCell(ICell);
          const Real DivCellExact = Setup.exactDivVec(X, Y);

          // Errors
          LInfCell(ICell)      = std::abs(DivCellNum - DivCellExact);
          LInfScaleCell(ICell) = std::abs(DivCellExact);
          L2Cell(ICell) = AreaCell(ICell) * LInfCell(ICell) * LInfCell(ICell);
          L2ScaleCell(ICell) =
              AreaCell(ICell) * LInfScaleCell(ICell) * LInfScaleCell(ICell);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfCell);
   const Real L2ErrorLoc   = sum(L2Cell);
   const Real LInfScaleLoc = maxVal(LInfScaleCell);
   const Real L2ScaleLoc   = sum(L2ScaleCell);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedDivErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedDivErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

int testGradient(Real RTol) {
   int Err;
   TestSetup Setup;

   const auto &mesh = HorzMesh::getDefault();
#ifdef HORZOPERATORS_TEST_PLANE
   const auto XCell = createDeviceMirrorCopy(mesh->XCellH);
   const auto YCell = createDeviceMirrorCopy(mesh->YCellH);
#else
   const auto XCell = createDeviceMirrorCopy(mesh->LonCellH);
   const auto YCell = createDeviceMirrorCopy(mesh->LatCellH);
#endif

   // Prepare operator input
   Array1DReal ScalarCell("ScalarCell", mesh->NCellsSize);
   parallelFor(
       {mesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          const Real X      = XCell(ICell);
          const Real Y      = YCell(ICell);
          ScalarCell(ICell) = Setup.exactScalar(X, Y);
       });

   // Perform halo exchange
   Halo MyHalo(MachEnv::getDefaultEnv(), Decomp::getDefault());
   auto ScalarCellH = createHostMirrorCopy(ScalarCell);
   MyHalo.exchangeFullArrayHalo(ScalarCellH, OnCell);
   deepCopy(ScalarCell, ScalarCellH);

#ifdef HORZOPERATORS_TEST_PLANE
   const auto XEdge = createDeviceMirrorCopy(mesh->XEdgeH);
   const auto YEdge = createDeviceMirrorCopy(mesh->YEdgeH);
#else
   const auto XEdge = createDeviceMirrorCopy(mesh->LonEdgeH);
   const auto YEdge = createDeviceMirrorCopy(mesh->LatEdgeH);
#endif
   const auto &AngleEdge = mesh->AngleEdge;
   const auto &DcEdge    = mesh->DcEdge;
   const auto &DvEdge    = mesh->DvEdge;

   // Compute element-wise errors
   Array1DReal LInfEdge("LInfEdge", mesh->NEdgesOwned);
   Array1DReal L2Edge("L2Edge", mesh->NEdgesOwned);
   Array1DReal LInfScaleEdge("LInfScaleEdge", mesh->NEdgesOwned);
   Array1DReal L2ScaleEdge("L2ScaleEdge", mesh->NEdgesOwned);
   GradientOnEdge GradientEdge(mesh);
   parallelFor(
       {mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          // Numerical result
          const Real GradScalarNum = GradientEdge(IEdge, ScalarCell);

          // Exact result
          const Real X                = XEdge(IEdge);
          const Real Y                = YEdge(IEdge);
          const Real GradScalarExactX = Setup.exactGradScalarX(X, Y);
          const Real GradScalarExactY = Setup.exactGradScalarY(X, Y);
          const Real EdgeNormalX      = std::cos(AngleEdge(IEdge));
          const Real EdgeNormalY      = std::sin(AngleEdge(IEdge));
          const Real GradScalarExact =
              EdgeNormalX * GradScalarExactX + EdgeNormalY * GradScalarExactY;

          LInfEdge(IEdge)      = std::abs(GradScalarNum - GradScalarExact);
          LInfScaleEdge(IEdge) = std::abs(GradScalarExact);
          const Real AreaEdge  = DcEdge(IEdge) * DvEdge(IEdge) / 2;
          L2Edge(IEdge)        = AreaEdge * LInfEdge(IEdge) * LInfEdge(IEdge);
          L2ScaleEdge(IEdge) =
              AreaEdge * LInfScaleEdge(IEdge) * LInfScaleEdge(IEdge);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfEdge);
   const Real LInfScaleLoc = maxVal(LInfScaleEdge);
   const Real L2ErrorLoc   = sum(L2Edge);
   const Real L2ScaleLoc   = sum(L2ScaleEdge);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedGradErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedGradErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

int testCurl(Real RTol) {
   int Err;
   TestSetup Setup;
   const auto &mesh = HorzMesh::getDefault();

#ifdef HORZOPERATORS_TEST_PLANE
   const auto XEdge = createDeviceMirrorCopy(mesh->XEdgeH);
   const auto YEdge = createDeviceMirrorCopy(mesh->YEdgeH);
#else
   const auto XEdge = createDeviceMirrorCopy(mesh->LonEdgeH);
   const auto YEdge = createDeviceMirrorCopy(mesh->LatEdgeH);
#endif
   const auto &AngleEdge = mesh->AngleEdge;

   // Prepare operator input
   Array1DReal VecEdge("VecEdge", mesh->NEdgesSize);
   parallelFor(
       {mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          const Real X = XEdge(IEdge);
          const Real Y = YEdge(IEdge);

          const Real VecExactX   = Setup.exactVecX(X, Y);
          const Real VecExactY   = Setup.exactVecY(X, Y);
          const Real EdgeNormalX = std::cos(AngleEdge(IEdge));
          const Real EdgeNormalY = std::sin(AngleEdge(IEdge));
          VecEdge(IEdge) = EdgeNormalX * VecExactX + EdgeNormalY * VecExactY;
       });

   // Perform halo exchange
   Halo MyHalo(MachEnv::getDefaultEnv(), Decomp::getDefault());
   auto VecEdgeH = createHostMirrorCopy(VecEdge);
   MyHalo.exchangeFullArrayHalo(VecEdgeH, OnEdge);
   deepCopy(VecEdge, VecEdgeH);

#ifdef HORZOPERATORS_TEST_PLANE
   const auto XVertex = createDeviceMirrorCopy(mesh->XVertexH);
   const auto YVertex = createDeviceMirrorCopy(mesh->YVertexH);
#else
   const auto XVertex = createDeviceMirrorCopy(mesh->LonVertexH);
   const auto YVertex = createDeviceMirrorCopy(mesh->LatVertexH);
#endif
   const auto &AreaTriangle = mesh->AreaTriangle;

   // Compute element-wise errors
   Array1DReal LInfVertex("LInfVertex", mesh->NVerticesOwned);
   Array1DReal LInfScaleVertex("LInfScaleVertex", mesh->NVerticesOwned);
   Array1DReal L2Vertex("L2Vertex", mesh->NVerticesOwned);
   Array1DReal L2ScaleVertex("L2ScaleVertex", mesh->NVerticesOwned);
   CurlOnVertex CurlVertex(mesh);
   parallelFor(
       {mesh->NVerticesOwned}, KOKKOS_LAMBDA(int IVertex) {
          // Numerical result
          const Real CurlNum = CurlVertex(IVertex, VecEdge);

          // Exact result
          const Real X         = XVertex(IVertex);
          const Real Y         = YVertex(IVertex);
          const Real CurlExact = Setup.exactCurlVec(X, Y);

          // Errors
          LInfVertex(IVertex)      = std::abs(CurlNum - CurlExact);
          LInfScaleVertex(IVertex) = std::abs(CurlExact);
          L2Vertex(IVertex) =
              AreaTriangle(IVertex) * LInfVertex(IVertex) * LInfVertex(IVertex);
          L2ScaleVertex(IVertex) = AreaTriangle(IVertex) *
                                   LInfScaleVertex(IVertex) *
                                   LInfScaleVertex(IVertex);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfVertex);
   const Real LInfScaleLoc = maxVal(LInfScaleVertex);
   const Real L2ErrorLoc   = sum(L2Vertex);
   const Real L2ScaleLoc   = sum(L2ScaleVertex);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedCurlErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedCurlErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

int testRecon(Real RTol) {
   int Err;
   TestSetup Setup;

   const auto &mesh = HorzMesh::getDefault();
#ifdef HORZOPERATORS_TEST_PLANE
   const auto XEdge = createDeviceMirrorCopy(mesh->XEdgeH);
   const auto YEdge = createDeviceMirrorCopy(mesh->YEdgeH);
#else
   const auto XEdge = createDeviceMirrorCopy(mesh->LonEdgeH);
   const auto YEdge = createDeviceMirrorCopy(mesh->LatEdgeH);
#endif
   const auto &AngleEdge = mesh->AngleEdge;

   // Prepare operator input
   Array1DReal VecEdge("VecEdge", mesh->NEdgesSize);
   parallelFor(
       {mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          const Real X = XEdge(IEdge);
          const Real Y = YEdge(IEdge);

          const Real VecExactX   = Setup.exactVecX(X, Y);
          const Real VecExactY   = Setup.exactVecY(X, Y);
          const Real EdgeNormalX = std::cos(AngleEdge(IEdge));
          const Real EdgeNormalY = std::sin(AngleEdge(IEdge));
          VecEdge(IEdge) = EdgeNormalX * VecExactX + EdgeNormalY * VecExactY;
       });

   // Perform halo exchange
   Halo MyHalo(MachEnv::getDefaultEnv(), Decomp::getDefault());
   auto VecEdgeH = createHostMirrorCopy(VecEdge);
   MyHalo.exchangeFullArrayHalo(VecEdgeH, OnEdge);
   deepCopy(VecEdge, VecEdgeH);

   const auto &DcEdge = mesh->DcEdge;
   const auto &DvEdge = mesh->DvEdge;

   // Compute element-wise errors
   Array1DReal LInfEdge("LInfEdge", mesh->NEdgesOwned);
   Array1DReal LInfScaleEdge("LInfScaleEdge", mesh->NEdgesOwned);
   Array1DReal L2Edge("L2Edge", mesh->NEdgesOwned);
   Array1DReal L2ScaleEdge("L2ScaleEdge", mesh->NEdgesOwned);
   TangentialReconOnEdge TanReconEdge(mesh);
   parallelFor(
       {mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          // Numerical result
          const Real VecReconNum = TanReconEdge(IEdge, VecEdge);

          // Exact result
          const Real X             = XEdge(IEdge);
          const Real Y             = YEdge(IEdge);
          const Real VecX          = Setup.exactVecX(X, Y);
          const Real VecY          = Setup.exactVecY(X, Y);
          const Real EdgeTangentX  = -std::sin(AngleEdge(IEdge));
          const Real EdgeTangentY  = std::cos(AngleEdge(IEdge));
          const Real VecReconExact = EdgeTangentX * VecX + EdgeTangentY * VecY;

          // Errors
          LInfEdge(IEdge)      = std::abs(VecReconNum - VecReconExact);
          LInfScaleEdge(IEdge) = std::abs(VecReconExact);
          const Real AreaEdge  = DcEdge(IEdge) * DvEdge(IEdge) / 2;
          L2Edge(IEdge)        = AreaEdge * LInfEdge(IEdge) * LInfEdge(IEdge);
          L2ScaleEdge(IEdge) =
              AreaEdge * LInfScaleEdge(IEdge) * LInfScaleEdge(IEdge);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfEdge);
   const Real LInfScaleLoc = maxVal(LInfScaleEdge);
   const Real L2ErrorLoc   = sum(L2Edge);
   const Real L2ScaleLoc   = sum(L2ScaleEdge);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedReconErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedReconErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

//------------------------------------------------------------------------------
// The initialization routine for Operators testing
int initOperatorsTest(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefaultEnv();
   MPI_Comm DefComm = DefEnv->getComm();

   Err = IO::init(DefComm);
   if (Err != 0) {
      LOG_ERROR("HorzMeshTest: error initializing parallel IO");
   }

   Err = Decomp::init();
   if (Err != 0) {
      LOG_ERROR("HorzMeshTest: error initializing default decomposition");
   }

   Err = HorzMesh::init();
   if (Err != 0) {
      LOG_ERROR("HorzMeshTest: error initializing default mesh");
   }

   return Err;
}

void finalizeOperatorsTest() {
   HorzMesh::clear();
   Decomp::clear();
   MachEnv::removeAll();
   Kokkos::finalize();
   MPI_Finalize();
}

int main(int argc, char *argv[]) {
   int Err = initOperatorsTest(argc, argv);
   if (Err != 0) {
      LOG_CRITICAL("OperatorsTest: Error initializing");
   }

   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 1e-10;

   int DivErr = testDivergence(RTol);
   if (DivErr == 0) {
      LOG_INFO("OperatorsTest: Divergence PASS");
   } else {
      Err = DivErr;
      LOG_INFO("OperatorsTest: Divergence FAIL");
   }

   int GradErr = testGradient(RTol);
   if (GradErr == 0) {
      LOG_INFO("OperatorsTest: Gradient PASS");
   } else {
      Err = GradErr;
      LOG_INFO("OperatorsTest: Gradient FAIL");
   }

   int CurlErr = testCurl(RTol);
   if (CurlErr == 0) {
      LOG_INFO("OperatorsTest: Curl PASS");
   } else {
      Err = CurlErr;
      LOG_INFO("OperatorsTest: Curl FAIL");
   }

   int ReconErr = testRecon(RTol);
   if (Err == 0) {
      LOG_INFO("OperatorsTest: Recon PASS");
   } else {
      Err = ReconErr;
      LOG_INFO("OperatorsTest: Recon FAIL");
   }

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Successful completion");
   }

   finalizeOperatorsTest();
} // end of main
//===-----------------------------------------------------------------------===/