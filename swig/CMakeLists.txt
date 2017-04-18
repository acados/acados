cmake_minimum_required(VERSION 2.8.11)

find_package(SWIG REQUIRED)
# SWIG wrapper file introduces some warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers -Wno-deprecated-register")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers")
endif ()

include(${SWIG_USE_FILE})

# Python
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

find_package(Casadi)

# SWIG does not support target_include_directories()
include_directories("${PROJECT_SOURCE_DIR}")
set_source_files_properties(acados.i PROPERTIES CPLUSPLUS ON )
# SWIG typemaps debug flag
# set(CMAKE_SWIG_FLAGS "-debug-tmused")
swig_add_module(acados python acados.i)
swig_link_libraries(acados ${PYTHON_LIBRARIES} acados qpoases blasfeo qpdunes)
if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    swig_link_libraries(acados ooqpgensparse ooqpsparse ooqpgondzio ooqpbase ma27 blas ${FLIB} m)
endif()

install(TARGETS ${SWIG_MODULE_acados_REAL_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${PROJECT_BINARY_DIR}/swig/acados.py DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES list_of_arrays.py DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES test_ocp_nlp.py DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES test_ocp_qp.py DESTINATION ${CMAKE_INSTALL_PREFIX})

add_test(NAME py_test_ocp_nlp COMMAND "${CMAKE_COMMAND}" -E chdir ${CMAKE_INSTALL_PREFIX} python3 test_ocp_nlp.py)
add_test(NAME py_test_ocp_qp COMMAND "${CMAKE_COMMAND}" -E chdir ${CMAKE_INSTALL_PREFIX} python3 test_ocp_qp.py)
