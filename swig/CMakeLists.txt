cmake_minimum_required(VERSION 2.8.11)

find_package(SWIG REQUIRED)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
execute_process(COMMAND python-config --prefix OUTPUT_VARIABLE python_prefix)
string(STRIP ${python_prefix} python_prefix)
set(PYTHON_LIBRARIES "${python_prefix}/lib/libpython2.7.dylib")

include(${SWIG_USE_FILE})

# Python
find_package(NumPy REQUIRED)
include_directories("/usr/local/include/python2.7")
include_directories("/usr/local/lib/python2.7/site-packages/numpy/core/include")
make_directory(${PROJECT_BINARY_DIR}/python/acados)
set_source_files_properties(acados.i  PROPERTIES  CPLUSPLUS ON)
set(CMAKE_SWIG_FLAGS "-debug-tmsearch")
swig_add_module(acados python acados.i)
swig_link_libraries(acados ${PYTHON_LIBRARIES} acados)
