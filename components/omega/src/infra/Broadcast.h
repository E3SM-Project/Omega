#ifndef OMEGA_BROADCAST_H
#define OMEGA_BROADCAST_H
//===-- infra/Broadcast.h - Broadcasting values --*- C++ -*-===//
//
/// \file
/// \brief Defines MPI broadcasting functions
///
/// This header defines functions to broadcast values from one MPI rank
/// to another.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"
#include "mpi.h"

namespace OMEGA {

// blocking broadcast scalar
int Broadcast(I4 &value, const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast = 0);
int Broadcast(I4 &value, const int rankBcast);

int Broadcast(I8 &value, const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast = 0);
int Broadcast(I8 &value, const int rankBcast);

int Broadcast(R4 &value, const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast = 0);
int Broadcast(R4 &value, const int rankBcast);

int Broadcast(R8 &value, const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast = 0);
int Broadcast(R8 &value, const int rankBcast);

int Broadcast(bool &value, const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast = 0);
int Broadcast(bool &value, const int rankBcast);

int Broadcast(std::string &value,
              const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast  = 0);
int Broadcast(std::string &value, const int rankBcast);

// blocking broadcast array
int Broadcast(std::vector<I4> &value,
              const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast  = 0);
int Broadcast(std::vector<I4> &value, const int rankBcast);

int Broadcast(std::vector<I8> &value,
              const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast  = 0);
int Broadcast(std::vector<I8> &value, const int rankBcast);

int Broadcast(std::vector<R4> &value,
              const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast  = 0);
int Broadcast(std::vector<R4> &value, const int rankBcast);

int Broadcast(std::vector<R8> &value,
              const MachEnv *InEnv = MachEnv::getDefaultEnv(),
              const int rankBcast  = 0);
int Broadcast(std::vector<R8> &value, const int rankBcast);

// NOTE: Elements of vector<bool> seem to be non-addressable
// int Broadcast(std::vector<bool> &value,
//              const MachEnv *InEnv = MachEnv::getDefaultEnv(),
//              const int rankBcast = 0);
// int Broadcast(std::vector<bool> &value, const int rankBcast);

/* To be implemented

// non-blocking broadcast scalar
void IBroadcast(const I4 value, const MachEnv *InEnv, const int rankIbcast = 0);
void IBroadcast(const I8 value, const MachEnv *InEnv, const int rankIbcast = 0);
void IBroadcast(const R4 value, const MachEnv *InEnv, const int rankIbcast = 0);
void IBroadcast(const R8 value, const MachEnv *InEnv, const int rankIbcast = 0);
void IBroadcast(const Real value, const MachEnv *InEnv, const int rankIbcast =
0); void IBroadcast(const bool value, const MachEnv *InEnv, const int rankIbcast
= 0); void IBroadcast(const std::string value, const MachEnv *InEnv, const int
rankIbcast = 0);

// non-blocking broadcast array
void IBroadcast(const std::vector<I4> value, const MachEnv *InEnv, const int
rankBcast = 0); void IBroadcast(const std::vector<I8> value, const MachEnv
*InEnv, const int rankBcast = 0); void IBroadcast(const std::vector<R4> value,
const MachEnv *InEnv, const int rankBcast = 0); void IBroadcast(const
std::vector<R8> value, const MachEnv *InEnv, const int rankBcast = 0); void
IBroadcast(const std::vector<Real> value, const MachEnv *InEnv, const int
rankBcast = 0); void IBroadcast(const std::vector<bool> value, const MachEnv
*InEnv, const int rankBcast = 0);

*/

} // namespace OMEGA

#endif // OMEGA_BROADCAST_H
