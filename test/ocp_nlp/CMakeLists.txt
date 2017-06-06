set(TEST_OCP_NLP_GENERATED_FILES
    "${PROJECT_BINARY_DIR}/test/ocp_nlp/jac_chain_nm2.c"
    "${PROJECT_BINARY_DIR}/test/ocp_nlp/jac_chain_nm3.c"
    "${PROJECT_BINARY_DIR}/test/ocp_nlp/jac_chain_nm4.c"
    "${PROJECT_BINARY_DIR}/test/ocp_nlp/vde_chain_nm2.c"
    "${PROJECT_BINARY_DIR}/test/ocp_nlp/vde_chain_nm3.c"
    "${PROJECT_BINARY_DIR}/test/ocp_nlp/vde_chain_nm4.c"
)

set(GENERATION_DIR "${CMAKE_SOURCE_DIR}/test/ocp_nlp/chain/casadi")
set(SED_COMMAND "'s/real_t\ a0/mem\ =\ 0;\ mem\ +=\ 0;\ w\ =\ 0;\ w\ +=\ 0;\ iw\ =\ 0;\ iw\ +=\ 0;  real_t\ a0/g'")

add_custom_command(OUTPUT ${TEST_OCP_NLP_GENERATED_FILES}
    COMMAND "octave-cli" "--no-gui" "--path" "${GENERATION_DIR}\:${EXTERNAL_SRC_DIR}/casadi-octave-v3.1.1"
        "${GENERATION_DIR}/generate_chain_mass.m"
    COMMAND "sed" "-i.bak" "-e" "${SED_COMMAND}" ${TEST_OCP_NLP_GENERATED_FILES}
    COMMENT "Generating data for OCP NLP test"
)

add_library(ocp_nlp_gen OBJECT ${TEST_OCP_NLP_GENERATED_FILES})
