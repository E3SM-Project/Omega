# build Omega

# Add source files for the library
file(GLOB LIBSRC_FILES ocn/*.cpp)
#set(LIBSRC_FILES ocn/miniWeather.cpp)

# Create the library target
add_library(${OMEGA_LIB_NAME} ${LIBSRC_FILES})
set_target_properties(${OMEGA_LIB_NAME} PROPERTIES COMPILE_FLAGS "${OMEGA_DEV_MINIWEATHER_EXE_DEFS}")

# Link the spdlog and YAKL libraries
target_link_libraries(${OMEGA_LIB_NAME} spdlog "${OMEGA_LDFLAGS}")
target_include_directories(${OMEGA_LIB_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${OMEGA_LIB_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/infra")

# include yakl cmake utility
include(${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/YAKL/yakl_utils.cmake)

# Enable caching for incremental compilation
set_property(TARGET ${OMEGA_LIB_NAME} PROPERTY CMAKE_CACHE_FILE_DEPENDS "CMakeCache.txt")

# help the library target to be built with yakl
yakl_process_target(${OMEGA_LIB_NAME})

# handles cuda cases
if (YAKL_ARCH STREQUAL "CUDA")
  set_target_properties(${OMEGA_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)
  if (CMAKE_VERSION VERSION_GREATER "3.18.0")
    set_target_properties(${OMEGA_LIB_NAME} PROPERTIES CUDA_ARCHITECTURES OFF)
  endif()
endif()

# build Omega executable
if(OMEGA_BUILD_EXECUTABLE)
  set(EXESRC_FILES
    drivers/main.cpp
  )
  # Create the executable target
  add_executable(${OMEGA_EXE_NAME} ${EXESRC_FILES})
  set_target_properties(${OMEGA_EXE_NAME} PROPERTIES COMPILE_FLAGS "${OMEGA_DEV_MINIWEATHER_EXE_DEFS}")

  # Link the Omega library
  target_link_libraries(${OMEGA_EXE_NAME} ${OMEGA_LIB_NAME})
  target_include_directories(${OMEGA_EXE_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_include_directories(${OMEGA_EXE_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/infra")

  # help the executable target to be built with yakl
  yakl_process_target(${OMEGA_EXE_NAME})

  # handles cuda cases
  if (YAKL_ARCH STREQUAL "CUDA")
    set_target_properties(${OMEGA_EXE_NAME} PROPERTIES LINKER_LANGUAGE CXX)
    if (CMAKE_VERSION VERSION_GREATER "3.18.0")
      set_target_properties(${OMEGA_EXE_NAME} PROPERTIES CUDA_ARCHITECTURES OFF)
    endif()
  endif()

  # Enable caching for incremental compilation
  set_property(TARGET ${OMEGA_EXE_NAME} PROPERTY CMAKE_CACHE_FILE_DEPENDS "CMakeCache.txt")

endif()
