include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    collections
    GIT_REPOSITORY "https://github.com/jchristopherson/collections"
    GIT_TAG main
)

FetchContent_MakeAvailable(collections)

set_target_properties(collections PROPERTIES
    OUTPUT_NAME "collections"
    NO_SONAME ON
    LIBRARY_OUTPUT_DIRECTORY ${MT_LIBS_DIR}
)

set(collections_INCLUDE_DIR ${collections_BINARY_DIR}/include)
set(collections_INCLUDE_DIR ${collections_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(collections_LIBRARY collections)
set(collections_LIBRARY ${collections_LIBRARY} PARENT_SCOPE)