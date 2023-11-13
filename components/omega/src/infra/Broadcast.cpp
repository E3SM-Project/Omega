//===-- infra/Broacast.cpp - Implements Broadcasting functions --*- C++ -*-===//
//
/// \file
/// \brief implements blocking and non-blocking broadcasting functions
///
/// This implements blocking and non-blocking broadcasting functions. Two
/// function names are overloaded accross Omega data types: I4, I8, R4, R8,
/// Real, bool, and std::string: 1) Broadcast for blocking mode and 2) mpiIbcast
/// for non- blocking mode.
//
//===----------------------------------------------------------------------===//

#include "Broadcast.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Broadcast I4 scalar value
int Broadcast(I4 &value, const MachEnv *InEnv, const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast(&value, 1, MPI_INT32_T, rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast for I4 data type

int Broadcast(I4 &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast for I4 data type

//------------------------------------------------------------------------------
// Broadcast I8 scalar value
int Broadcast(I8 &value, const MachEnv *InEnv, const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast(&value, 1, MPI_INT64_T, rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast for I8 data type

int Broadcast(I8 &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast for I8 data type

//------------------------------------------------------------------------------
// Broadcast R4 scalar value
int Broadcast(R4 &value, const MachEnv *InEnv, const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast(&value, 1, MPI_FLOAT, rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(R4 &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R8 scalar value
int Broadcast(R8 &value, const MachEnv *InEnv, const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast(&value, 1, MPI_DOUBLE, rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(R8 &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast bool scalar value
int Broadcast(bool &value, const MachEnv *InEnv, const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast(&value, 1, MPI_C_BOOL, rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(bool &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast std::string value
int Broadcast(std::string &value, const MachEnv *InEnv, const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast((void *)value.c_str(), value.size(), MPI_CHAR, rankBcast,
                      InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(std::string &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast I4 array
int Broadcast(std::vector<I4> &value, const MachEnv *InEnv,
              const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast((void *)value.data(), value.size(), MPI_INT32_T,
                      rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(std::vector<I4> &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast I8 array
int Broadcast(std::vector<I8> &value, const MachEnv *InEnv,
              const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast((void *)value.data(), value.size(), MPI_INT64_T,
                      rankBcast, InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(std::vector<I8> &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R4 array
int Broadcast(std::vector<R4> &value, const MachEnv *InEnv,
              const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast((void *)value.data(), value.size(), MPI_FLOAT, rankBcast,
                      InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(std::vector<R4> &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R8 array
int Broadcast(std::vector<R8> &value, const MachEnv *InEnv,
              const int rankBcast) {
   int retval = 0;

   retval = MPI_Bcast((void *)value.data(), value.size(), MPI_DOUBLE, rankBcast,
                      InEnv->getComm());

   return retval;
} // end Broadcast

int Broadcast(std::vector<R8> &value, const int rankBcast) {
   return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
}

//------------------------------------------------------------------------------
// Broadcast bool array
// Elements of vector<bool> seem to be non-addressable
// int Broadcast(std::vector<bool> &value, const MachEnv *InEnv, const int
// rankBcast ) {
//    int retval = 0;
//
//    retval = MPI_Bcast((void *)value.data(), value.size(), MPI_C_BOOL,
//    rankBcast, InEnv->getComm());
//
//    return retval;
//} // end Broadcast
//
// int Broadcast(std::vector<bool> &value, const int rankBcast
//) {
//    return Broadcast(value, MachEnv::getDefaultEnv(), rankBcast);
//}

} // namespace OMEGA
