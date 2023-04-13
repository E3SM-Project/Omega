# Add external packages

# Add the spdlog header-only library
add_library(spdlog INTERFACE)
target_include_directories(spdlog INTERFACE spdlog/include)

# Add the YAKL library
add_subdirectory(YAKL)
if (CMAKE_VERSION VERSION_GREATER "3.18.0")
  set_property(TARGET yakl PROPERTY CUDA_ARCHITECTURES OFF)
endif()
include_directories(YAKL)
