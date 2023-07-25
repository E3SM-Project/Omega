//===-- base/MachEnv.cpp - machine environment methods ----------*- C++ -*-===//
//
// The machine environment defines a number of parameters associated with
// the message-passing and threading environments. It also can describe
// the node-level hardware environment, including any useful accelerator
// or processor-level constants. Multiple machine environments can be defined,
// mostly associated with subsets of the processor decomposition, but
// a default machine env must be created very early in model
// initialization. All environments can be retrieved by name, though a
// specific retrieval for the most common default environment is provided
// for efficiency.
//
//===----------------------------------------------------------------------===//

#include "MachEnv.h"
#include "mpi.h"

#include <string>
#include <map>
// Note that we should replace iostream and std::cerr with the logging
// capability once that is enabled.
#include <iostream>

namespace OMEGA {

// create the static class members
MachEnv* MachEnv::DefaultEnv = nullptr;
std::map<std::string, MachEnv*> MachEnv::AllEnvs;

// Constructors
//------------------------------------------------------------------------------
// Default constructor that fills a MachEnv with mostly invalid values.
// Environments are instead created using the relevant create function.

MachEnv::MachEnv(){

   // Set the communicator to the MPI equivalent of a null ptr
   Comm = MPI_COMM_NULL;

   // Set all other values to invalid numbers
   MyTask         = -999;
   NumTasks       = -999;
   MasterTask     = -999;
   MasterTaskFlag = false;
   MemberFlag     = true;
   NumThreads     = -999;

} // end default constructor

//------------------------------------------------------------------------------
// Add a MachEnv to list of instances
//
void MachEnv::addEnv(const std::string Name, // [in] name of environment
                     MachEnv* NewEnv         // [in] environment to add
                     ){

   // Check to see if an environment of the same name already exists and
   // if so, exit with an error
   if (AllEnvs.find(Name) != AllEnvs.end()){
      std::cerr << "Attempted to create a MachEnv with name " << Name <<
                   " but an Env of that name already exists ";
      return;
   }

   // Add to list of environments
   AllEnvs[Name] = NewEnv;

} // end addEnv

//------------------------------------------------------------------------------
// Initializes the Machine Environment by creating the DefaultEnv for Omega

void MachEnv::init(const MPI_Comm InComm // [in] communicator to use
                  ){

   // Create the Env structure and make it persistent (static)
   static MachEnv DefEnv;

   // Set the communicator to the input communicator by duplicating it
   MPI_Comm_dup(InComm, &(DefEnv.Comm));

   // get task ID and set local MPI task
   MPI_Comm_rank(DefEnv.Comm, &(DefEnv.MyTask));

   // get total number of MPI tasks
   MPI_Comm_size(DefEnv.Comm, &(DefEnv.NumTasks));

   // Set task 0 as master
   DefEnv.MasterTask = 0;

   // determine if this task is the master task
   if (DefEnv.MyTask == DefEnv.MasterTask){
      DefEnv.MasterTaskFlag = true;
   } else {
      DefEnv.MasterTaskFlag = false;
   }

   // All tasks are members of this communicator's group
   DefEnv.MemberFlag = true;

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   DefEnv.NumThreads = omp_get_num_threads();
#else
   DefEnv.NumThreads = 1;
#endif

   // Add to list of environments
   OMEGA::MachEnv::addEnv("Default", &DefEnv);

   // Retrieve this environment and set pointer to DefaultEnv
   MachEnv::DefaultEnv = getEnv("Default");

} // end init Mach Env

// Create functions
//------------------------------------------------------------------------------
// Create a new environment from a contiguous subset of an
// existing environment starting from same root task.

void MachEnv::createEnv(const std::string Name, // [in] name of environment
                        const MachEnv* InEnv,   // [in] parent MachEnv
                        const int NewSize       // [in] num tasks in new env
                        ){

   // Create the new Env structure and make it persistent (static)
   static MachEnv NewEnv;

   // Get communicator from old environment
   MPI_Comm InComm = InEnv->getComm();

   // Error checks on new size
   int OldSize;
   MPI_Comm_size(InComm, &OldSize);
   if (NewSize > OldSize){
      std::cerr << "Invalid NewSize in MachEnv constructor " <<
                   "NewSize = " << NewSize << 
                   " is larger than old size" << OldSize << std::endl;
      return;
   }

   // First retrieve the group associated with the old environment
   MPI_Group InGroup;
   MPI_Comm_group(InComm, &InGroup);

   // Define the range of tasks in new group (0, NewSize-1)
   int NRanges = 1;
   int LastTask = NewSize - 1;
   int Range[1][3] = {0, LastTask, 1};

   // Create a new group with the new range
   MPI_Group NewGroup;
   MPI_Group_range_incl(InGroup, NRanges, Range, &NewGroup);

   // Create the communicator for the new group
   MPI_Comm_create(InComm, NewGroup, &NewEnv.Comm);

   // determine whether this task is a part of the new communicator
   if (NewEnv.Comm != MPI_COMM_NULL){

      NewEnv.MemberFlag = true;
      // get task ID for local MPI task
      MPI_Comm_rank(NewEnv.Comm, &NewEnv.MyTask);
      // get total number of MPI tasks
      MPI_Comm_size(NewEnv.Comm, &NewEnv.NumTasks);

   } else {

      NewEnv.MemberFlag = false;
      NewEnv.MyTask     = -999;
      NewEnv.NumTasks   = -999;
   }

   // Set task 0 as master
   NewEnv.MasterTask = 0;

   // determine if this task is the master task
   if (NewEnv.MyTask == NewEnv.MasterTask){
      NewEnv.MasterTaskFlag = true;
   } else {
      NewEnv.MasterTaskFlag = false;
   }

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NewEnv.NumThreads = omp_get_num_threads();
#else
   NewEnv.NumThreads = 1;
#endif

   // Add to list of environments
   OMEGA::MachEnv::addEnv(Name, &NewEnv);

} // end create with contiguous range

//------------------------------------------------------------------------------
// Create a new environment from a strided subset of an
// existing environment

void MachEnv::createEnv(const std::string Name, // [in] name of environment
                        const MachEnv* InEnv,   // [in] parent MachEnv
                        const int NewSize,      // [in] num tasks in new env
                        const int Begin,        // [in] starting parent task
                        const int Stride        // [in] stride for tasks to incl
                        ){

   // Create the new Env structure and make it persistent (static)
   static MachEnv NewEnv;

   // First retrieve the group associated with the parent environment
   MPI_Comm InComm = InEnv->getComm();
   MPI_Group InGroup;
   MPI_Comm_group(InComm, &InGroup);

   // Define the range of tasks in new group based on the input
   // start task and stride
   int NRanges = 1;
   int LastTask = Begin + (NewSize-1)*Stride;
   int Range[1][3] = {Begin, LastTask, Stride};

   // Error checks for valid range
   int OldSize;
   MPI_Comm_size(InComm, &OldSize);
   if (LastTask >= OldSize){
      std::cerr << "Invalid range in strided MachEnv constructor " <<
                   "LastTask = " << LastTask << 
                   " is larger than old size" << OldSize << std::endl;
      return;
   }

   // Create a new group with the new range
   MPI_Group NewGroup;
   MPI_Group_range_incl(InGroup, NRanges, Range, &NewGroup);

   // Create the communicator for the new group
   MPI_Comm_create(InComm, NewGroup, &(NewEnv.Comm));

   // determine whether this task is a part of the new communicator
   // and set quantities accordingly
   if (NewEnv.Comm != MPI_COMM_NULL){ // this task is part of the group

      NewEnv.MemberFlag = true;
      // get task ID for local MPI task
      MPI_Comm_rank(NewEnv.Comm, &(NewEnv.MyTask));
      // get total number of MPI tasks
      MPI_Comm_size(NewEnv.Comm, &(NewEnv.NumTasks));

   } else {
      NewEnv.MemberFlag = false;
      NewEnv.MyTask     = -999;
      NewEnv.NumTasks   = -999;
   }

   // Set task 0 as master
   NewEnv.MasterTask = 0;

   // determine if this task is the master task
   if (NewEnv.MyTask == NewEnv.MasterTask){
      NewEnv.MasterTaskFlag = true;
   } else {
      NewEnv.MasterTaskFlag = false;
   }

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NewEnv.NumThreads = omp_get_num_threads();
#else
   NewEnv.NumThreads = 1;
#endif

   // Add to list of environments
   OMEGA::MachEnv::addEnv(Name, &NewEnv);

} // end constructor using strided range

//------------------------------------------------------------------------------
// Create a new environment from a custom subset of an
// existing environment, supplying list of parent tasks to include

void MachEnv::createEnv(
                 const std::string Name, // [in] name of environment
                 const MachEnv* InEnv,   // [in] parent MPI communicator
                 const int NewSize,      // [in] num tasks in new env
                 const int Tasks[]       // [in] vector of parent tasks to incl
                 ){

   // Create the new Env structure and make it persistent (static)
   static MachEnv NewEnv;

   // Error checks on valid tasks
   MPI_Comm InComm = InEnv->getComm();
   int OldSize;
   int ThisTask;
   MPI_Comm_size(InComm, &OldSize);
   for (int i=0; i < NewSize; ++i){
      ThisTask = Tasks[i];
      if (ThisTask < 0 || ThisTask >= OldSize){
         std::cerr << "Invalid task in MachEnv constructor " <<
                      "task = " << ThisTask << 
                      " is < 0 or larger than old size" << OldSize
                      << std::endl;
      return;
      }
   }

   // First retrieve the group associated with the input communicator
   MPI_Group InGroup;
   MPI_Comm_group(InComm, &InGroup);

   // Create a new group with the selected tasks
   MPI_Group NewGroup;
   MPI_Group_incl(InGroup, NewSize, Tasks, &NewGroup);

   // Create the communicator for the new group
   MPI_Comm_create(InComm, NewGroup, &(NewEnv.Comm));

   // determine whether this task is a part of the new communicator
   // and set values accordingly
   if (NewEnv.Comm != MPI_COMM_NULL){ // member of new group

      NewEnv.MemberFlag = true;
      // get task ID for local MPI task
      MPI_Comm_rank(NewEnv.Comm, &(NewEnv.MyTask));
      // get total number of MPI tasks
      MPI_Comm_size(NewEnv.Comm, &(NewEnv.NumTasks));

   } else {
      NewEnv.MemberFlag = false;
      NewEnv.MyTask     = -999;
      NewEnv.NumTasks   = -999;
   }

   // Set task 0 as master
   NewEnv.MasterTask = 0;

   // determine if this task is the master task
   if (NewEnv.MyTask == NewEnv.MasterTask){
      NewEnv.MasterTaskFlag = true;
   } else {
      NewEnv.MasterTaskFlag = false;
   }

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NewEnv.NumThreads = omp_get_num_threads();
#else
   NewEnv.NumThreads = 1;
#endif

   // Add to list of environments
   OMEGA::MachEnv::addEnv(Name, &NewEnv);


} // end constructor with selected tasks

// Remove/delete functions
//------------------------------------------------------------------------------
// Remove environment

void MachEnv::removeEnv (const std::string name // [in] name of env to remove
      ){

   AllEnvs.erase(name);

} // end removeEnv

// Retrieval functions
//------------------------------------------------------------------------------
// Get default environment
MachEnv* MachEnv::getDefaultEnv() { return MachEnv::DefaultEnv; }

//------------------------------------------------------------------------------
// Get environment by name
MachEnv* MachEnv::getEnv(
                const std::string name ///< [in] name of environment
                ){

    // look for an instance of this name
    auto it = AllEnvs.find(name);

    // if found, return the environment pointer
    if (it != AllEnvs.end()){
       return it->second;

    // otherwise print an error and return a null pointer
    } else {
       std::cerr << "Attempt to retrieve a non-existent MachEnv " <<
                    name << " has not been defined or has been removed";
       return nullptr;
    }

} // end getEnv

//------------------------------------------------------------------------------
// Get communicator for an environment
MPI_Comm MachEnv::getComm() const { return Comm; }

//------------------------------------------------------------------------------
// Get local task/rank ID
int MachEnv::getMyTask() const { return MyTask; }

//------------------------------------------------------------------------------
// Get total number of MPI tasks/ranks
int MachEnv::getNumTasks() const { return NumTasks; }

//------------------------------------------------------------------------------
// Get task ID for the master task (typically 0)
int MachEnv::getMasterTask() const { return MasterTask; }

//------------------------------------------------------------------------------
// Determine whether local task is the master

bool MachEnv::isMasterTask() const { return MasterTaskFlag; }

//------------------------------------------------------------------------------
// Determine whether local task is in this communicator's group

bool MachEnv::isMember() const { return MemberFlag; }

//------------------------------------------------------------------------------
// Set task ID for the master task (if not 0)

int MachEnv::setMasterTask(const int TaskID){

   int Err=0;

   if (TaskID >= 0 && TaskID < NumTasks){
      MasterTask = TaskID;
      if (MyTask == MasterTask){
         MasterTaskFlag = true;
      } else {
         MasterTaskFlag = false;
      }
   } else {
      std::cerr << "Error: invalid TaskID sent to MachEnv.setMasterTask" <<
                   std::endl;
      Err = -1;
   }
   return Err;

} // end setMasterTask

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
