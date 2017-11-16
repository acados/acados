find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
find_package(NumPy REQUIRED)
find_package(CasadiPython REQUIRED)

# do not install casadi twice
if(CASADIPYTHON_FOUND)
    add_library(casadi_python SHARED IMPORTED)
    set_property(TARGET casadi_python PROPERTY IMPORTED_LOCATION "${CASADI_PYTHON_LIBRARY}")
    include_directories("${CASADI_PYTHON_INCLUDE_DIR}")
    get_filename_component(CASADI_PYTHON_LIBRARY_DIR "${CASADI_PYTHON_LIBRARY}" DIRECTORY)
    list(APPEND CMAKE_INSTALL_RPATH "${CASADI_PYTHON_LIBRARY_DIR}")
endif()

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

set_source_files_properties(../acados.i PROPERTIES CPLUSPLUS ON)
if(ACADOS_WITH_OOQP)
    set_source_files_properties(../acados.i PROPERTIES SWIG_FLAGS -DOOQP)
endif()

swig_add_module(acados python ../acados.i)
swig_link_libraries(acados ${PYTHON_LIBRARIES} acados casadi_python)

if(ACADOS_WITH_OOQP)
    target_compile_definitions(${SWIG_MODULE_acados_REAL_NAME} PRIVATE OOQP)
endif()

install(TARGETS ${SWIG_MODULE_acados_REAL_NAME} DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/acados.py DESTINATION lib)