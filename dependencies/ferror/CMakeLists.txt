# Get the macros and functions we'll need
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    ferror
    GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
)

FetchContent_MakeAvailable(ferror)

set(ferror_INCLUDE_DIR ${ferror_BINARY_DIR}/include)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)
configure_file(
    "${ferror_SOURCE_DIR}/include/ferror.h"
    "${ferror_INCLUDE_DIR}/ferror.h"
    COPYONLY
)

# Make a parent-scope variable for the library
set(ferror_LIBRARY ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
