# Add source files for the library
set(SRC_FILES
  infra/logging.cpp
)

# Create the library target
add_library(OmegaLib ${SRC_FILES})

# Link the spdlog library
target_link_libraries(OmegaLib spdlog)

#
## Create the library target
#add_library(MyLargeSoftwareLib ${SRC_FILES})
#
## Configure compiler-specific and machine-specific information
#if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
#  # Configure for GCC compiler
#  target_compile_options(MyLargeSoftwareLib PRIVATE "-Wall")
#elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "MSVC")
#  # Configure for MSVC compiler
#  target_compile_options(MyLargeSoftwareLib PRIVATE "/W4")
#endif()

# Enable caching for incremental compilation
set_property(TARGET OmegaLib PROPERTY CMAKE_CACHE_FILE_DEPENDS "CMakeCache.txt")

# Set output directory
set_target_properties(OmegaLib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${OMEGA_OUTPUT_DIRECTORY}")
set_target_properties(OmegaLib PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${OMEGA_OUTPUT_DIRECTORY}")

if(OMEGA_BUILD_EXECUTABLE)
  # Create the executable target
  add_executable(omega.exe drivers/main.cpp)

  # Link the spdlog library
  target_link_libraries(omega.exe OmegaLib)

  # Enable caching for incremental compilation
  set_property(TARGET omega.exe PROPERTY CMAKE_CACHE_FILE_DEPENDS "CMakeCache.txt")

  # Set output directory
  set_target_properties(omega.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${OMEGA_OUTPUT_DIRECTORY}")
endif()

