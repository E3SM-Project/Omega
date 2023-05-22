# Add external packages

# Add the spdlog header-only library
add_library(spdlog INTERFACE)
target_include_directories(spdlog INTERFACE spdlog/include)

# Add the YAKL library
#add_subdirectory(YAKL)
add_subdirectory(${CMAKE_SOURCE_DIR}/../../externals/YAKL ${CMAKE_CURRENT_BINARY_DIR}/YAKL)
if (CMAKE_VERSION VERSION_GREATER "3.18.0")
  set_property(TARGET yakl PROPERTY CUDA_ARCHITECTURES OFF)
endif()
include_directories(${CMAKE_SOURCE_DIR}/../../externals/YAKL)
#include_directories(YAKL)
