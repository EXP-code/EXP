# Checks for the existence of feenableexcept() and sets the flag
# HAVE_FE_ENABLE if the code snippet compiles.  The HAVE_FE_ENABLE
# preprocessor flag will then be set in config_exp.h if TRUE.
#
# MDW 04/07/23

set(_fe_snippet
  "
  #include <cfenv>
  #include <limits>
  #include <sstream>
  int main()
  {
    feenableexcept(FE_DIVBYZERO|FE_INVALID);
    std::ostringstream description;
    const double lower_bound = -std::numeric_limits<double>::max();
    description << lower_bound;
    return 0;
  }
  "
  )

# Check the test source
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("${_fe_snippet}" FE_ENABLE_FOUND)

# Parse the response flag and provide config info
if(FE_ENABLE_FOUND)
  set(HAVE_FE_ENABLE TRUE)
  message(STATUS "We have <feenableexcept>: will compile FPE handlers.")
else()
  message(STATUS "We do NOT have <feenableexcept>: not compling FPE handlers.")
endif()
