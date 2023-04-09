# Add source files for the library
set(LIBSRC_FILES
  infra/logging.cpp
  ocn/miniWeather.cpp
)

# Create the library target
add_library(OmegaLib ${LIBSRC_FILES})
set_target_properties(OmegaLib PROPERTIES COMPILE_FLAGS "${OMEGA_DEV_MINIWEATHER_EXE_DEFS}")

# Link the spdlog and YAKL libraries
target_link_libraries(OmegaLib spdlog "${LDFLAGS}")
target_include_directories(OmegaLib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

include(${CMAKE_CURRENT_SOURCE_DIR}/../external/YAKL/yakl_utils.cmake)

# Enable caching for incremental compilation
set_property(TARGET OmegaLib PROPERTY CMAKE_CACHE_FILE_DEPENDS "CMakeCache.txt")

yakl_process_target(OmegaLib)

if ("${YAKL_ARCH}" STREQUAL "CUDA")
  set_target_properties(OmegaLib PROPERTIES LINKER_LANGUAGE CXX)
  if (${CMAKE_VERSION} VERSION_GREATER "3.18.0")
    set_target_properties(OmegaLib PROPERTIES CUDA_ARCHITECTURES OFF)
  endif()
endif()

# Set output directory
set_target_properties(OmegaLib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${OMEGA_OUTPUT_DIRECTORY}")
set_target_properties(OmegaLib PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${OMEGA_OUTPUT_DIRECTORY}")

if(OMEGA_BUILD_EXECUTABLE)
  # Create the executable target
  add_executable(omega.exe drivers/main.cpp)
  set_target_properties(omega.exe PROPERTIES COMPILE_FLAGS "${OMEGA_DEV_MINIWEATHER_EXE_DEFS}")

  # Link the Omega library
  target_link_libraries(omega.exe OmegaLib)
  target_include_directories(omega.exe PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

  yakl_process_target(omega.exe)

  if ("${YAKL_ARCH}" STREQUAL "CUDA")
    set_target_properties(omega.exe PROPERTIES LINKER_LANGUAGE CXX)
    if (${CMAKE_VERSION} VERSION_GREATER "3.18.0")
      set_target_properties(omega.exe PROPERTIES CUDA_ARCHITECTURES OFF)
    endif()
  endif()

  # Enable caching for incremental compilation
  set_property(TARGET omega.exe PROPERTY CMAKE_CACHE_FILE_DEPENDS "CMakeCache.txt")

  # Set output directory
  set_target_properties(omega.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${OMEGA_OUTPUT_DIRECTORY}")
endif()

