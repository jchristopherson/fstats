include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fplot
    GIT_REPOSITORY "https://github.com/jchristopherson/fplot"
)

FetchContent_MakeAvailable(fplot)

set(fplot_INCLUDE_DIR ${fplot_BINARY_DIR}/include)
set(fplot_INCLUDE_DIR ${fplot_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fplot_LIBRARY fplot)
set(fplot_LIBRARY ${fplot_LIBRARY} PARENT_SCOPE)
