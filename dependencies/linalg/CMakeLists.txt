# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# If found, use the installed version; else, import the library
if (${linalg_FOUND})
    # Inform the user of what's going on
    message(STATUS "LINALG (${linalg_VERSION}) library found.")

    # Get the mod file location
    get_target_property(linalg_INCLUDE_DIR linalg INTERFACE_INCLUDE_DIRECTORIES)
else()
    # Inform the user of what's going on
    message(STATUS "LINALG library not found.  Downloading appropriate repository.")
    
    # Fetch the proper content
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
        OVERRIDE_FIND_PACKAGE
    )

    FetchContent_MakeAvailable(linalg)

    if (WIN32)
        if (BUILD_SHARED_LIBS)
            add_custom_command(
                TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:linalg>
                $<TARGET_FILE_DIR:${PROJECT_NAME}
            )
        endif()
    endif()
    
    set(linalg_INCLUDE_DIR ${linalg_BINARY_DIR}/include)
endif()

# Make a parent-scope variable for the library
set(linalg_LIBRARY linalg)
set(linalg_LIBRARY ${linalg_LIBRARY} PARENT_SCOPE)

# Make a parent-scope variable locating the include directory for linalg
set(linalg_INCLUDE_DIR ${linalg_INCLUDE_DIR} PARENT_SCOPE)