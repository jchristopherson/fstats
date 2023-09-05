include(FetchContent)
FetchContent_Declare(
    fortran-csv-module
    GIT_TAG "origin/master"
    GIT_REPOSITORY https://github.com/jchristopherson/fortran-csv-module
)
FetchContent_MakeAvailable(fortran-csv-module)
set(fortran-csv-module_INCLUDE_DIR ${fortran-csv-module_BINARY_DIR}/fcsv_include)
set(fortran-csv-module_LIBRARY fcsv)

set(fortran-csv-module_INCLUDE_DIR ${fortran-csv-module_INCLUDE_DIR} PARENT_SCOPE)
set(fortran-csv-module_LIBRARY ${fortran-csv-module_LIBRARY} PARENT_SCOPE)