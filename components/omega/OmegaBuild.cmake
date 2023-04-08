# define build control variables(BCVs)

# Omega build steps
# * -1: not started
# * 0: started
# * 100: finished setup step
# * 200: finished update and check step
# * 300: finished build step
# * 400: finished output step

####################################
# Create internal Omega build-contr#
# olling variables                 #
####################################

# Omega build step: -1: not started, 0:started, n: step number
set(OMEGA_BUILD_STEP "-1")


####################################
# Create user-controllable Omega bu#
# ild-controlling variables        #
####################################
macro(setup_common_variables)

  # Omega build output directory
  if(NOT OMEGA_OUTPUT_DIRECTORY)
    set(OMEGA_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bld")
  endif()


endmacro()

macro(setup_standalone_build)

  setup_common_variables()

  set(OMEGA_BUILD_STEP 100)
endmacro()

macro(setup_e3sm_build)

  setup_common_variables()

  set(OMEGA_BUILD_STEP 100)
endmacro()

macro(update_cmake_variables)

endmacro()

macro(check_setup)

  if(${OMEGA_BUILD_STEP} LESS 0 OR ${OMEGA_BUILD_STEP} GREATER 100)
    message(FATAL_ERROR "OMEGA_BUILD_STEP should be 1 after setup step, but is ${OMEGA_BUILD_STEP}.")
  endif()

  set(OMEGA_BUILD_STEP 200)

endmacro()

# fullfill the requrests from build initiators
macro(wrap_outputs)
endmacro()

