include(FetchContent)

# Get FERROR
add_subdirectory(ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)

# Get LINALG
add_subdirectory(linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)

# Get COLLECTIONS
add_subdirectory(collections)
set(collections_LIBRARY ${collections_LIBRARY} PARENT_SCOPE)
set(collections_INCLUDE_DIR ${collections_INCLUDE_DIR} PARENT_SCOPE)