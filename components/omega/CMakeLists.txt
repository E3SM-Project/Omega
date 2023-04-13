# Omega Top-level CMakeLists.txt
#
# - supports two build modes:
#    1) standalone build
#    2) e3sm build
#
# - variable name prefix: "OMEGA_"
#     * controls the Omega building


###########################################################
# STEP 1: Setup                                           #
#                                                         #
# Assign values to the Omega build-controlling variables  #
###########################################################

# define variables and macros for Omega build
include(${CMAKE_CURRENT_SOURCE_DIR}/OmegaBuild.cmake)

# handle build types 
if (NOT DEFINED PROJECT_NAME)
  # enter standa-lone build

  cmake_minimum_required(VERSION 3.16)
  project(${OMEGA_PROJECT_NAME})

  set(CMAKE_CXX_STANDARD 11)

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

###########################################################
# STEP 2: Update                                          #
#                                                         #
# update cmake & yakl variables, and adjust               #
# build-controlling variables                             #
###########################################################

# update CMake variables from the build-controlling variables
update_variables()

###########################################
# place holder for further build-controls #
###########################################

# check if all-required variables have proper values
check_setup()


###########################################################
# STEP 3: Build                                           #
#                                                         #
# build Omega, and generate tests                          #
###########################################################

add_subdirectory(external)
add_subdirectory(src)

if(OMEGA_BUILD_TEST)

  enable_testing()
  add_subdirectory(test)

endif()

###########################################################
# STEP 4: Output                                          #
#                                                         #
# prepare outputs according to the output-controlling     #
# variables                                               #
###########################################################

wrap_outputs()
