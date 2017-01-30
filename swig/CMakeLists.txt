cmake_minimum_required(VERSION 2.8.11)

find_package(SWIG REQUIRED)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
execute_process(COMMAND python-config --prefix OUTPUT_VARIABLE python_prefix)
string(STRIP ${python_prefix} python_prefix)
set(PYTHON_LIBRARIES "${python_prefix}/lib/libpython2.7.dylib")

include(${SWIG_USE_FILE})

# Python
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
find_package(NumPy 1.8 REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
make_directory(${PROJECT_BINARY_DIR}/python/acados)
set_source_files_properties(acados.i  PROPERTIES  CPLUSPLUS ON)
# set(CMAKE_SWIG_FLAGS "-debug-tmsearch")
swig_add_module(acados python acados.i)
swig_link_libraries(acados ${PYTHON_LIBRARIES} acados)
