
# Add a test for the stand-alone executable
add_test(NAME OMEGA_TEST COMMAND ${OMEGA_EXE_NAME})

# Configure test dependencies
set_tests_properties(OMEGA_TEST PROPERTIES DEPENDS ${OMEGA_EXE_NAME})
