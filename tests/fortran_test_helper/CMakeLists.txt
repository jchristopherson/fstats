include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fortran_test_helper
    GIT_REPOSITORY "https://github.com/jchristopherson/fortran_test_helper"
    GIT_TAG main
)

FetchContent_MakeAvailable(fortran_test_helper)

set(fortran_test_helper_INCLUDE_DIR ${fortran_test_helper_BINARY_DIR}/include)
set(fortran_test_helper_INCLUDE_DIR ${fortran_test_helper_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fortran_test_helper_LIBRARY fortran_test_helper)
set(fortran_test_helper_LIBRARY ${fortran_test_helper_LIBRARY} PARENT_SCOPE)