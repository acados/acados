find_package(Matlab REQUIRED)
find_package(CasadiMatlab REQUIRED)

if(CASADIMATLAB_FOUND)
    add_library(casadi_matlab SHARED IMPORTED)
    set(CASADI_LIBRARY_NAME "${CMAKE_SHARED_LIBRARY_PREFIX}casadi${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(CASADI_LIBRARY_FILE "${CASADI_MATLAB_LIBRARY}/${CASADI_LIBRARY_NAME}")
    set_property(TARGET casadi_matlab PROPERTY IMPORTED_LOCATION "${CASADI_LIBRARY_FILE}")
    include_directories("${CASADI_MATLAB_INCLUDE_DIR}")
    list(APPEND CMAKE_INSTALL_RPATH "${CASADI_MATLAB_LIBRARY}")
endif()

set_source_files_properties(../acados.i PROPERTIES CPLUSPLUS ON)

swig_module_initialize(acadosMEX matlab)
swig_add_source_to_module(acadosMEX FALSE MEX_FILE ../acados.i)

if(EXISTS ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
    add_library(acadosMATLAB MODULE ${MEX_FILE} ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
else()
    add_library(acadosMATLAB MODULE ${MEX_FILE})
endif()
target_include_directories(acadosMATLAB PUBLIC ${MATLAB_INCLUDE_DIR})
target_link_libraries(acadosMATLAB ${MATLAB_LIBRARIES} acados casadi_matlab)
set_target_properties(acadosMATLAB
    PROPERTIES
        OUTPUT_NAME "acadosMEX"
        PREFIX ""
        SUFFIX .${MATLAB_MEX_EXT}
)
install(TARGETS acadosMATLAB DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+acados DESTINATION lib)
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SwigGet.m
        ${CMAKE_CURRENT_BINARY_DIR}/SwigMem.m
        ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m
    DESTINATION lib
)
