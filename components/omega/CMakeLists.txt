# Omega Top-level CMakeLists.txt
#
# - supports two build modes:
#    1) standalone and
#    2) embedded(e3sm)
#
# - uses a prefix of "OMEGA_" in all
#    Omega build controlling variables(OBCVs)


####################################
# STEP 1: Setup                    #
#                                  #
# Assign proper values to build-con#
# trolling variables               #
####################################

# include variables and macros for Omega build
include(${CMAKE_CURRENT_SOURCE_DIR}/OmegaBuild.cmake)

# handle multiple types of Omega build
if (NOT DEFINED PROJECT_NAME)
  # standalone build

  cmake_minimum_required(VERSION 3.16)
  project(OmegaOceanModel)

  set(CMAKE_CXX_STANDARD 11)

  # Set the build type to Release by default
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
  endif()

  # update variables for standalone build
  setup_standalone_build()

else()
  # enter embedded build

  if(DEFINED E3SM_SOURCE_DIR)

    # update variables for e3sm build
    setup_e3sm_build()

  else()
    message(FATAL_ERROR "UNKNOWN OMEGA BUILD MODE")

  endif()

endif()

####################################
# STEP 2: Update                   #
#                                  #
# modify build-controlling variable#
# s to accomodate further changes  #
####################################

# update CMake variables from Omega build-controlling variables(OBCVs)
update_cmake_variables()

# check if all-required variables have proper values
check_setup()

####################################
# STEP 3: Build                    #
#                                  #
# build Omega and generate tests   #
####################################

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(test)

####################################
# STEP 4: Output                   #
#                                  #
# prepare outputs according to the #
# output-controlling variables     #
####################################

wrap_outputs()
