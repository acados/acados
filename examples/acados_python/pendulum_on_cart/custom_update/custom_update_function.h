/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#include "acados_solver_pendulum.h"


// Called at the end of solver creation.
// This is allowed to allocate memory and store the pointer to it into capsule->custom_update_memory.
int custom_update_init_function(pendulum_solver_capsule* capsule);


// Custom update function that can be called between solver calls
int custom_update_function(pendulum_solver_capsule* capsule, double* data, int data_len);

// Called just before destroying the solver.
// Responsible to free allocated memory, stored at capsule->custom_update_memory.
int custom_update_terminate_function(pendulum_solver_capsule* capsule);
