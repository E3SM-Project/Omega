###########################
# Internal variables      #
###########################

set(OMEGA_PROJECT_NAME "OmegaOceanModel")
set(OMEGA_EXE_NAME "omega.exe")
set(OMEGA_LIB_NAME "OmegaLib")


###########################
# Public variables        #
###########################
macro(setup_common_variables)

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

  # prints generates all cmake variables
  #get_cmake_property(_variableNames VARIABLES)
  #list (SORT _variableNames)
  #foreach (_variableName ${_variableNames})
  #    message(STATUS "${_variableName}=${${_variableName}}")
  #endforeach()

  setup_common_variables()

endmacro()


################################
# Set cmake and YAKL variables #
################################
macro(update_variables)

  if(OMEGA_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${OMEGA_INSTALL_PREFIX})
  endif()

  if(OMEGA_ARCH)
    set(YAKL_ARCH ${OMEGA_ARCH})
  
    if(OMEGA_${OMEGA_ARCH}_FLAGS)
      set(YAKL_${OMEGA_ARCH}_FLAGS ${OMEGA_${OMEGA_ARCH}_FLAGS})
    endif()

  endif()

endmacro()



################################
# Verify variable integrity    #
################################
macro(check_setup)

endmacro()



################################
# Prepare output               #
################################
macro(wrap_outputs)

  if(OMEGA_INSTALL_PREFIX)

    install(TARGETS ${OMEGA_LIB_NAME} LIBRARY DESTINATION "${OMEGA_INSTALL_PREFIX}/lib")
  
    if(OMEGA_BUILD_EXECUTABLE)
      install(TARGETS ${OMEGA_EXE_NAME} RUNTIME DESTINATION "${OMEGA_INSTALL_PREFIX}/bin")
    endif()

  endif()

endmacro()

