find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
find_package(NumPy REQUIRED)
find_package(Casadi REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

# set(CMAKE_SWIG_FLAGS "-keyword")

set_source_files_properties(../acados.i PROPERTIES CPLUSPLUS ON)

swig_add_module(acados python ../acados.i)
swig_link_libraries(acados ${PYTHON_LIBRARIES} acados casadi)

install(TARGETS ${SWIG_MODULE_acados_REAL_NAME} DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/acados.py DESTINATION lib)