# Check that exp can be run
add_test(NAME expRunTest
  COMMAND exp -v WORKING_DIRECTORY src)

# Check that pyEXP can be loaded
add_test(NAME pyexpLoadTest
  COMMAND python3 test.py
  WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")

# Append the pyEXP  build directory to PYTHONPATH
set_property(TEST pyexpLoadTest PROPERTY ENVIRONMENT_MODIFICATION
  "PYTHONPATH=path_list_prepend:${CMAKE_BINARY_DIR}/pyexp}")

# Look for a file
add_test(NAME scriptTest
  COMMAND "${CMAKE_CURRENT_LIST_DIR}/script_test.sh"
  WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")
