set(Matlab_ROOT_DIR "$ENV{MATLAB_ROOT}")
find_package(Matlab COMPONENTS MX_LIBRARY ENG_LIBRARY MAT_LIBRARY REQUIRED)
find_package(CasadiMatlab REQUIRED)

if(CASADIMATLAB_FOUND)
    add_library(casadi_matlab UNKNOWN IMPORTED)
    set_property(TARGET casadi_matlab PROPERTY IMPORTED_LOCATION "${CASADI_MATLAB_LIBRARY}")
    include_directories("${CASADI_MATLAB_INCLUDE_DIR}")
    get_filename_component(CASADI_MATLAB_LIBRARY_DIR "${CASADI_MATLAB_LIBRARY}" DIRECTORY)
    list(APPEND CMAKE_INSTALL_RPATH "${CASADI_MATLAB_LIBRARY_DIR}")
endif()

set_source_files_properties(../acados.i PROPERTIES CPLUSPLUS ON)
if(ACADOS_WITH_OOQP)
    set_source_files_properties(../acados.i PROPERTIES SWIG_FLAGS -DOOQP)
endif()

if(CMAKE_SYSTEM MATCHES "Windows")
    if(${CMAKE_GENERATOR} MATCHES "MinGW Makefiles")
        set_source_files_properties(../acados.i PROPERTIES SWIG_FLAGS "-DSWIG_WIN_MINGW")
    endif()
endif()

swig_module_initialize(acadosMEX matlab)
swig_add_source_to_module(acadosMEX FALSE MEX_FILE ../acados.i)

if(EXISTS ${Matlab_INCLUDE_DIRS}/../version/cpp_mexapi_version.cpp)
    add_library(acadosMATLAB MODULE ${MEX_FILE} ${Matlab_INCLUDE_DIRS}/../version/cpp_mexapi_version.cpp)
else()
    add_library(acadosMATLAB MODULE ${MEX_FILE})
endif()
target_include_directories(acadosMATLAB PUBLIC ${Matlab_INCLUDE_DIRS})
target_link_libraries(acadosMATLAB acados casadi_matlab ${Matlab_LIBRARIES})
set_target_properties(acadosMATLAB
    PROPERTIES
        OUTPUT_NAME "acadosMEX"
        PREFIX ""
        SUFFIX .${Matlab_MEX_EXTENSION}
)

if(ACADOS_WITH_OOQP)
    target_compile_definitions(acadosMATLAB PRIVATE OOQP)
endif()

install(TARGETS acadosMATLAB DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+acados DESTINATION lib)
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SwigGet.m
        ${CMAKE_CURRENT_BINARY_DIR}/SwigMem.m
        ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m
    DESTINATION lib
)
