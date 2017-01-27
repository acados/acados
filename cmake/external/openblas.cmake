if(NOT OpenBLAS_LIB)
    find_package(OpenBLAS)
endif()

add_library(openblas SHARED IMPORTED)
set_property(TARGET openblas PROPERTY IMPORTED_LOCATION "${OpenBLAS_LIB}")
