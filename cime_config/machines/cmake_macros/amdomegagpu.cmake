string(APPEND CMAKE_C_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CPPDEFS_DEBUG " -DYAKL_DEBUG")

if (compile_threaded)
  string(APPEND CMAKE_Fortran_FLAGS   " -fopenmp")
  string(APPEND CMAKE_C_FLAGS   " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS  " -fopenmp")
endif()

string(APPEND CPPDEFS " -DNO_R16 -DCPRAMD -DFORTRANUNDERSCORE")

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")
