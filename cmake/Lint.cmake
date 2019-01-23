# Linter
find_package(PythonInterp 3)

# XXX temporarely disable lint !!!!!!!!!!!
set(LINT_COMMAND echo lint is disabled!)
add_custom_target(lint ${LINT_COMMAND})
return()
# XXX end of temporarely disable lint !!!!!!!!!!!

if(NOT PYTHONINTERP_FOUND OR CMAKE_SYSTEM_NAME MATCHES "Windows" OR CMAKE_SYSTEM_NAME MATCHES "dSpace")
    message(STATUS "Disabling lint")
    return()
endif()

set(FIND_FILES_TO_LINT find acados examples interfaces swig test -type f -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" -o -name "*.i")
set(FIND_FILES_TO_LINT ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_FILES_TO_LINT})
execute_process(COMMAND ${FIND_FILES_TO_LINT} OUTPUT_VARIABLE FILES_TO_LINT)
string(REPLACE "\n" " " FILES_TO_LINT ${FILES_TO_LINT})
separate_arguments(FILES_TO_LINT)
set(LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --quiet --counting=detailed --extensions=c,cpp,h,hpp,i ${FILES_TO_LINT})
add_custom_target(lint ${LINT_COMMAND})
