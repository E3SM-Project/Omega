# define build control variables(BCVs)


##############################################################
# Create internal Omega build-controlling variables          #
##############################################################


##############################################################
# Create user-controllable Omega build-controlling variables #
##############################################################

macro(setup_common_variables)

  # build an executable target
  if(NOT OMEGA_BUILD_EXECUTABLE)
    set(OMEGA_BUILD_EXECUTABLE OFF)
  endif()

  # set output directory
  if(NOT OMEGA_OUTPUT_DIRECTORY)
    set(OMEGA_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bld")
  endif()


  # miniWeather macros for development
  if ("${OMEGA_DEV_MINIWEATHER_NX}" STREQUAL "")
    SET(OMEGA_DEV_MINIWEATHER_NX 100)
  endif()
  if ("${OMEGA_DEV_MINIWEATHER_NZ}" STREQUAL "")
    SET(OMEGA_DEV_MINIWEATHER_NZ 50)
  endif()
  if ("${OMEGA_DEV_MINIWEATHER_SIM_TIME}" STREQUAL "")
    SET(OMEGA_DEV_MINIWEATHER_SIM_TIME 1000)
  endif()
  if ("${OMEGA_DEV_MINIWEATHER_OUT_FREQ}" STREQUAL "")
    SET(OMEGA_DEV_MINIWEATHER_OUT_FREQ 10)
  endif()
  if ("${OMEGA_DEV_MINIWEATHER_DATA_SPEC}" STREQUAL "")
    SET(OMEGA_DEV_MINIWEATHER_DATA_SPEC DATA_SPEC_THERMAL)
  endif()

  SET(OMEGA_DEV_MINIWEATHER_EXE_DEFS "-D_NX=${OMEGA_DEV_MINIWEATHER_NX} -D_NZ=${OMEGA_DEV_MINIWEATHER_NZ} -D_SIM_TIME=${OMEGA_DEV_MINIWEATHER_SIM_TIME} -D_OUT_FREQ=${OMEGA_DEV_MINIWEATHER_OUT_FREQ} -D_DATA_SPEC=${OMEGA_DEV_MINIWEATHER_DATA_SPEC}")
  
endmacro()

macro(setup_standalone_build)

  setup_common_variables()

endmacro()

macro(setup_e3sm_build)

  setup_common_variables()

endmacro()

macro(update_cmake_variables)

endmacro()

macro(check_setup)

endmacro()

# fullfill the requrests from build initiators
macro(wrap_outputs)
endmacro()

