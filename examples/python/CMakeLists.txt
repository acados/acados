add_test(NAME ocp_nlp_py
    COMMAND ${PYTHON_EXECUTABLE} ${EXAMPLES_SRC_DIR}/python/ocp_nlp.py)

add_test(NAME ocp_qp_py
    COMMAND ${PYTHON_EXECUTABLE} ${EXAMPLES_SRC_DIR}/python/ocp_qp.py)