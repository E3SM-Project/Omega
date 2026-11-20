#ifndef OMEGA_BINNEDSUMS_H
#define OMEGA_BINNEDSUMS_H
//===-- analysis/BinnedSums.h - reproducible binned sums --------*- C++ -*-===//
//
/// \file
/// \brief Defines BinnedSums, reproducible sums of mesh entities by bin
///
/// BinnedSums sums values over the owned mesh entities in each bin, at each
/// of a number of levels, so that the result is the same from run to run and
/// does not depend on how the mesh is divided among MPI tasks.
///
/// Sums accumulated with atomic adds depend on the order in which the values
/// arrive, which on a GPU changes from one run to the next. BinnedSums
/// instead lists each bin's entities once, in local index order, and splits
/// each list into segments of at most SegmentSize entities. A team per
/// segment and level reduces its entities in parallel, and a team per bin
/// and level then reduces its segment partials. Each team writes only its own
/// result, and a team reduction with a fixed team size (OMEGA_TEAMSIZE)
/// combines its threads' partial sums in the same order on every run. Both
/// stages keep the high and low parts of each double-double sum. Each task's
/// partial sums go to globalSumDD, which adds them across tasks in
/// double-double as well.
///
/// The double-double sums carry about 106 bits, so a different division of
/// the mesh among tasks changes the result only in the rare case that the
/// exact sum lies almost exactly halfway between two doubles. This is the
/// same reproducibility as globalSum. The order of the team reductions
/// depends on the Kokkos backend and team size, so results are not
/// guaranteed to match bit for bit across different hardware.
///
/// A caller builds the lists from each entity's bin, then computes the sums
/// with a device functor giving the value of an entity at a level:
/// \code
///   Sums.build(NumBins, BinOfEntity); // BinOfEntity is a HostArray1DI4
///   auto Partials = Sums.localSums(
///       NLevels, KOKKOS_LAMBDA(int Entity, int Bin, int Level) {...});
///   std::vector<R8> Global = globalSumDD(Partials, Comm);
/// \endcode
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Error.h"
#include "OmegaKokkos.h"
#include "Reductions.h"

#include <complex>
#include <vector>

namespace OMEGA {

/// Reproducible sums of entity values by bin and level
class BinnedSums {
 public:
   /// Maximum number of entities a team reduces before its partial sum is
   /// combined with those of the other segments of the same bin
   static constexpr I4 SegmentSize = 256;

   /// Builds the lists of entities in each bin. Entity I belongs to bin
   /// BinOfEntity[I]; entities with a bin outside [0, NumBins) belong to none.
   /// Each list holds its entities in increasing index order.
   void build(I4 InNumBins,                    ///< [in] number of bins
              const HostArray1DI4 &BinOfEntity ///< [in] bin of each entity
   ) {
      NumBins            = InNumBins;
      const I4 NEntities = static_cast<I4>(BinOfEntity.extent(0));

      // Count the entities in each bin
      HostArray1DI4 BinCount("BinnedSumsBinCount", NumBins);
      for (I4 I = 0; I < NEntities; ++I) {
         const I4 Bin = BinOfEntity(I);
         if (Bin >= 0 && Bin < NumBins)
            ++BinCount(Bin);
      }

      // Offsets of each bin's list, and of each bin's first segment
      HostArray1DI4 BinStart("BinnedSumsBinStart", NumBins + 1);
      HostArray1DI4 BinSegStartH("BinnedSumsBinSegStart", NumBins + 1);
      for (I4 Bin = 0; Bin < NumBins; ++Bin) {
         BinStart(Bin + 1) = BinStart(Bin) + BinCount(Bin);
         const I4 NBinSegs = (BinCount(Bin) + SegmentSize - 1) / SegmentSize;
         BinSegStartH(Bin + 1) = BinSegStartH(Bin) + NBinSegs;
      }
      NSegments = BinSegStartH(NumBins);

      // Fill each bin's list in increasing entity order
      HostArray1DI4 EntitiesH("BinnedSumsEntities", BinStart(NumBins));
      HostArray1DI4 Next("BinnedSumsNext", NumBins);
      for (I4 Bin = 0; Bin < NumBins; ++Bin)
         Next(Bin) = BinStart(Bin);
      for (I4 I = 0; I < NEntities; ++I) {
         const I4 Bin = BinOfEntity(I);
         if (Bin >= 0 && Bin < NumBins)
            EntitiesH(Next(Bin)++) = I;
      }

      // Split each bin's list into segments
      HostArray1DI4 SegStartH("BinnedSumsSegStart", NSegments + 1);
      HostArray1DI4 SegBinH("BinnedSumsSegBin", NSegments);
      for (I4 Bin = 0; Bin < NumBins; ++Bin) {
         for (I4 Seg = BinSegStartH(Bin); Seg < BinSegStartH(Bin + 1); ++Seg) {
            SegStartH(Seg) =
                BinStart(Bin) + (Seg - BinSegStartH(Bin)) * SegmentSize;
            SegBinH(Seg) = Bin;
         }
      }
      SegStartH(NSegments) = BinStart(NumBins);

      Entities    = createDeviceMirrorCopy(EntitiesH);
      SegStart    = createDeviceMirrorCopy(SegStartH);
      SegBin      = createDeviceMirrorCopy(SegBinH);
      BinSegStart = createDeviceMirrorCopy(BinSegStartH);
      Built       = true;

      // Scratch arrays are sized on the next call to localSums
      ScratchLevels = -1;
   }

   /// Whether the lists have been built
   bool isBuilt() const { return Built; }

   /// Computes this task's double-double sums of Value(Entity, Bin, Level)
   /// over the entities in each bin, for each level in [0, NLevels). The
   /// result has NumBins * NLevels elements, Bin * NLevels + Level, each with
   /// the sum in its real part and the residual in its imaginary part, ready
   /// for globalSumDD.
   template <class F>
   std::vector<std::complex<double>>
   localSums(I4 NLevels,    ///< [in] number of levels
             const F &Value ///< [in] device functor (Entity, Bin, Level) -> R8
   ) {
      OMEGA_REQUIRE(Built, "BinnedSums: localSums called before build");

      if (ScratchLevels != NLevels) {
         ScratchLevels = NLevels;
         SegHi         = Array2DR8("BinnedSumsSegHi", NSegments, NLevels);
         SegLo         = Array2DR8("BinnedSumsSegLo", NSegments, NLevels);
         BinHi         = Array2DR8("BinnedSumsBinHi", NumBins, NLevels);
         BinLo         = Array2DR8("BinnedSumsBinLo", NumBins, NLevels);
      }

      OMEGA_SCOPE(LocEntities, Entities);
      OMEGA_SCOPE(LocSegStart, SegStart);
      OMEGA_SCOPE(LocSegBin, SegBin);
      OMEGA_SCOPE(LocBinSegStart, BinSegStart);
      OMEGA_SCOPE(LocSegHi, SegHi);
      OMEGA_SCOPE(LocSegLo, SegLo);
      OMEGA_SCOPE(LocBinHi, BinHi);
      OMEGA_SCOPE(LocBinLo, BinLo);
      const F LocValue = Value;

      parallelForOuter(
          "BinnedSumsSegments", {NSegments, NLevels},
          KOKKOS_LAMBDA(int Seg, int Level, const TeamMember &Team) {
             const I4 Bin = LocSegBin(Seg);
             DDValue Sum;
             parallelReduceInner(
                 Team, Range{LocSegStart(Seg), LocSegStart(Seg + 1) - 1},
                 INNER_LAMBDA(int J, DDValue &Accum) {
                    sumDDLocal(Accum.Hi, Accum.Lo,
                               LocValue(LocEntities(J), Bin, Level));
                 },
                 DDReducer(Sum));
             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    LocSegHi(Seg, Level) = Sum.Hi;
                    LocSegLo(Seg, Level) = Sum.Lo;
                 });
          });

      parallelForOuter(
          "BinnedSumsBins", {NumBins, NLevels},
          KOKKOS_LAMBDA(int Bin, int Level, const TeamMember &Team) {
             DDValue Sum;
             parallelReduceInner(
                 Team, Range{LocBinSegStart(Bin), LocBinSegStart(Bin + 1) - 1},
                 INNER_LAMBDA(int Seg, DDValue &Accum) {
                    sumDDPair(Accum.Hi, Accum.Lo, LocSegHi(Seg, Level),
                              LocSegLo(Seg, Level));
                 },
                 DDReducer(Sum));
             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    LocBinHi(Bin, Level) = Sum.Hi;
                    LocBinLo(Bin, Level) = Sum.Lo;
                 });
          });

      auto BinHiH = createHostMirrorCopy(BinHi);
      auto BinLoH = createHostMirrorCopy(BinLo);
      std::vector<std::complex<double>> Partials(NumBins * NLevels);
      for (I4 Bin = 0; Bin < NumBins; ++Bin)
         for (I4 Level = 0; Level < NLevels; ++Level)
            Partials[Bin * NLevels + Level] =
                std::complex<double>(BinHiH(Bin, Level), BinLoH(Bin, Level));
      return Partials;
   }

 private:
   struct DDValue {
      R8 Hi;
      R8 Lo;
   };

   /// Keeps the residual when team threads merge their partial sums.
   struct DDReducer {
      using reducer    = DDReducer;
      using value_type = DDValue;
      using result_view_type =
          Kokkos::View<value_type, HostMemSpace, Kokkos::MemoryUnmanaged>;

      KOKKOS_INLINE_FUNCTION
      explicit DDReducer(value_type &Value) : Result(&Value) {}

      KOKKOS_INLINE_FUNCTION void init(value_type &Value) const {
         Value.Hi = 0.0;
         Value.Lo = 0.0;
      }

      KOKKOS_INLINE_FUNCTION void join(value_type &Dest,
                                       const value_type &Src) const {
         sumDDPair(Dest.Hi, Dest.Lo, Src.Hi, Src.Lo);
      }

      KOKKOS_INLINE_FUNCTION value_type &reference() const {
         return *Result.data();
      }

      KOKKOS_INLINE_FUNCTION result_view_type view() const { return Result; }

      KOKKOS_INLINE_FUNCTION bool references_scalar() const { return true; }

      result_view_type Result;
   };

   bool Built       = false;
   I4 NumBins       = 0;
   I4 NSegments     = 0;
   I4 ScratchLevels = -1; ///< number of levels the scratch arrays hold
   Array1DI4 Entities;    ///< entities of each bin, in bin then index order
   Array1DI4 SegStart;    ///< offset in Entities of each segment (NSegments+1)
   Array1DI4 SegBin;      ///< bin of each segment
   Array1DI4 BinSegStart; ///< first segment of each bin (NumBins+1)
   Array2DR8 SegHi;       ///< segment sums, high parts
   Array2DR8 SegLo;       ///< segment sums, low parts
   Array2DR8 BinHi;       ///< bin sums, high parts
   Array2DR8 BinLo;       ///< bin sums, low parts
};

} // namespace OMEGA

#endif
