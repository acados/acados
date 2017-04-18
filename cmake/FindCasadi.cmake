execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import casadi; print(casadi.CasadiMeta_getVersion())"
    OUTPUT_VARIABLE CASADI_VERSION
)
string(STRIP "${CASADI_VERSION}" CASADI_VERSION)
if (NOT CASADI_VERSION)
    message(FATAL_ERROR "Casadi not found!")
endif ()
string(SUBSTRING "${CASADI_VERSION}" 0 1 CASADI_MAJOR_VERSION)
string(COMPARE EQUAL "${CASADI_MAJOR_VERSION}" "3" FOUND_CASADI_3)
if (NOT FOUND_CASADI_3)
    message(FATAL_ERROR "Casadi version 3 required. Found version: ${CASADI_VERSION}")
else ()
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c "import casadi, os, inspect; print(os.path.dirname(inspect.getfile(casadi)))"
        OUTPUT_VARIABLE Casadi_LIBRARY
        RESULT_VARIABLE Casadi_OUTPUT
    )
    string(STRIP "${Casadi_LIBRARY}" Casadi_LIBRARY)
    get_filename_component(Casadi_BASE_DIR "${Casadi_LIBRARY}/../" ABSOLUTE)
    find_path(Casadi_INCLUDE_DIR NAMES casadi/casadi.hpp
        HINTS ${Casadi_BASE_DIR}
        PATHS ${Casadi_BASE_DIR}/include
    )
    find_package_handle_standard_args(Casadi DEFAULT_MSG Casadi_LIBRARY Casadi_INCLUDE_DIR)
endif ()

