###########################
# Internal variables      #
###########################

set(OMEGA_PROJECT_NAME            "OmegaOceanModel")
set(OMEGA_EXE_NAME                "omega.exe")
set(OMEGA_LIB_NAME                "OmegaLib")

set(OMEGA_BUILD_MODES             "E3SM" "STANDALONE" "NOT_DEFINED")
set(OMEGA_BUILD_MODE              NOT_DEFINED CACHE STRING "Omega build mode")
set_property(CACHE OMEGA_BUILD_MODE PROPERTY STRINGS ${OMEGA_BUILD_MODES})

set(OMEGA_SOURCE_DIR              ${CMAKE_CURRENT_LIST_DIR})
set(OMEGA_DEFAULT_BUILD_TYPE      Release) # Debug or Release


###########################
# Public variables        #
###########################
macro(setup_common_variables)

  if(NOT DEFINED OMEGA_BUILD_TYPE)
    set(OMEGA_BUILD_TYPE ${OMEGA_DEFAULT_BUILD_TYPE})
  endif()

  if(NOT DEFINED OMEGA_E3SM_SOURCE_DIR)
    set(OMEGA_E3SM_SOURCE_DIR ${OMEGA_SOURCE_DIR}/..)
  endif()

  if(NOT DEFINED OMEGA_BUILD_TYPE)
    set(OMEGA_BUILD_TYPE ${OMEGA_DEFAULT_BUILD_TYPE})
  endif()

  # miniWeather macros for development
  if (NOT DEFINED OMEGA_DEV_MINIWEATHER_NX)
    SET(OMEGA_DEV_MINIWEATHER_NX 100)
  endif()
  if (NOT DEFINED OMEGA_DEV_MINIWEATHER_NZ)
    SET(OMEGA_DEV_MINIWEATHER_NZ 50)
  endif()
  if (NOT DEFINED OMEGA_DEV_MINIWEATHER_SIM_TIME)
    SET(OMEGA_DEV_MINIWEATHER_SIM_TIME 1000)
  endif()
  if (NOT DEFINED OMEGA_DEV_MINIWEATHER_OUT_FREQ)
    SET(OMEGA_DEV_MINIWEATHER_OUT_FREQ 10)
  endif()
  if (NOT DEFINED OMEGA_DEV_MINIWEATHER_DATA_SPEC)
    SET(OMEGA_DEV_MINIWEATHER_DATA_SPEC "DATA_SPEC_THERMAL")
  endif()

  SET(OMEGA_DEV_MINIWEATHER_EXE_DEFS "\
    -D_NX=${OMEGA_DEV_MINIWEATHER_NX} \
    -D_NZ=${OMEGA_DEV_MINIWEATHER_NZ} \
    -D_SIM_TIME=${OMEGA_DEV_MINIWEATHER_SIM_TIME} \
    -D_OUT_FREQ=${OMEGA_DEV_MINIWEATHER_OUT_FREQ} \
    -D_DATA_SPEC=${OMEGA_DEV_MINIWEATHER_DATA_SPEC}")
  
endmacro()


macro(setup_standalone_build)

  setup_common_variables()

  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/../omega AND
      EXISTS ${CMAKE_CURRENT_LIST_DIR}/../../components)

    set(E3SM_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../../components)

  else()
    # so far, we assume that Omega exists inside of E3SM.
    # However, we leave this else part for later usage.

  endif()

  set(OMEGA_BUILD_MODE "STANDALONE")

endmacro()

macro(setup_e3sm_build)

  # prints generates all cmake variables
  #get_cmake_property(_variableNames VARIABLES)
  #list (SORT _variableNames)
  #foreach (_variableName ${_variableNames})
  #    message(STATUS "${_variableName}=${${_variableName}}")
  #endforeach()

  set(OMEGA_E3SM_BINARY_DIR ${E3SM_BINARY_DIR})
  set(OMEGA_BUILD_TYPE ${E3SM_DEFAULT_BUILD_TYPE})


  setup_common_variables()

  #TODO: set OMEGA_ARCH according to E3SM variables
  set(OMEGA_ARCH "NOT_DEFINED")
  set(OMEGA_BUILD_MODE "E3SM")

endmacro()


################################
# Set cmake and YAKL variables #
################################
macro(update_variables)

  # Set the build type
  set(CMAKE_BUILD_TYPE ${OMEGA_BUILD_TYPE})

  if(OMEGA_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${OMEGA_INSTALL_PREFIX})
  endif()

  if(DEFINED OMEGA_ARCH)

    if(OMEGA_ARCH STREQUAL "NOT_DEFINED")
      set(YAKL_ARCH "")

    else()
      set(YAKL_ARCH ${OMEGA_ARCH})
 
      if(OMEGA_${YAKL_ARCH}_FLAGS)
        set(YAKL_${YAKL_ARCH}_FLAGS ${OMEGA_${YAKL_ARCH}_FLAGS})
      endif()

    endif()
 
  endif()

endmacro()



################################
# Verify variable integrity    #
################################
macro(check_setup)

  message("OMEGA_BUILD_MODE = ${OMEGA_BUILD_MODE}")

  if(OMEGA_BUILD_MODE STREQUAL "E3SM")
    message("In E3SM")

  elseif(${OMEGA_BUILD_MODE} STREQUAL "STANDALONE")
    message("In STANDALONE")

  else()

    message(FATAL_ERROR "OMEGA_BUILD_MODE is neither E3SM nor STANDALONE.")

  endif()

  if (NOT DEFINED YAKL_ARCH)
    message(FATAL_ERROR "YAKL_ARCH is not defined.")
  endif()

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

