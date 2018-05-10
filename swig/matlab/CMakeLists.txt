find_package(Matlab COMPONENTS MX_LIBRARY ENG_LIBRARY REQUIRED)
find_package(Casadi REQUIRED)

set_source_files_properties(../acados.i PROPERTIES CPLUSPLUS ON)

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
target_link_libraries(acadosMATLAB acados casadi ${Matlab_LIBRARIES})
set_target_properties(acadosMATLAB
    PROPERTIES
        OUTPUT_NAME "acadosMEX"
        PREFIX ""
        SUFFIX .${Matlab_MEX_EXTENSION}
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
