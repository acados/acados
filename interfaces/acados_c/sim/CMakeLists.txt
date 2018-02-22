
set(INTERFACES_ACADOS_C_SIM_SRC
${CMAKE_CURRENT_SOURCE_DIR}/sim_erk_integrator_i.c
${CMAKE_CURRENT_SOURCE_DIR}/sim_lifted_irk_integrator_i.c)

target_sources(acados PRIVATE ${INTERFACES_ACADOS_C_SIM_SRC})
