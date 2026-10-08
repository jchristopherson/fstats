include(FetchContent)
set(BUILD_TESTING OFF)
set(CMAKE_PROJECT_fortran-csv-module_INCLUDE "${CMAKE_CURRENT_LIST_DIR}/disable_csv_tests.cmake")
FetchContent_Declare(
    fortran-csv-module
    GIT_TAG "origin/master"
    GIT_REPOSITORY https://github.com/jchristopherson/fortran-csv-module
)
FetchContent_MakeAvailable(fortran-csv-module)
unset(CMAKE_PROJECT_fortran-csv-module_INCLUDE)
set_property(DIRECTORY "${fortran-csv-module_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL TRUE)
set(fortran-csv-module_INCLUDE_DIR ${fortran-csv-module_BINARY_DIR}/fcsv_include)
set(fortran-csv-module_LIBRARY fcsv)

set(fortran-csv-module_INCLUDE_DIR ${fortran-csv-module_INCLUDE_DIR} PARENT_SCOPE)
set(fortran-csv-module_LIBRARY ${fortran-csv-module_LIBRARY} PARENT_SCOPE)