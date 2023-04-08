## Add test source files
#set(TEST_FILES
#  test_file1.cpp
#  test_file2.cpp
#  ...
#)
#
## Create the test executable
#add_executable(MyTestExecutable ${TEST_FILES})
#
## Link with the library
#target_link_libraries(MyTestExecutable MyLargeSoftwareLib)
#
## Add tests using CTest
#enable_testing()
#
## Add a test for the stand-alone executable
#add_test(NAME StandAloneTest COMMAND MyLargeSoftware arg1 arg2)
#
## Add a test for the library mode
#add_test(NAME LibraryTest COMMAND MyTestExecutable)
#
## Configure test dependencies
#set_tests_properties(LibraryTest PROPERTIES DEPENDS MyTestExecutable)
#
