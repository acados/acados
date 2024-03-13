/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

// external
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include <blasfeo/include/blasfeo.h>

// acados
#include <acados/ocp_nlp/ocp_nlp_reg_common.h>
#include <acados/ocp_nlp/ocp_nlp_reg_mirror.h>
#include <acados/ocp_nlp/ocp_nlp_reg_project.h>
#include <acados/ocp_nlp/ocp_nlp_reg_convexify.h>
#include <acados/ocp_nlp/ocp_nlp_reg_noreg.h>



int main()
{
	printf("\nregularization example\n\n");

	int ii, jj;

    /************************************************
     * config
     ************************************************/

	acados_size_t config_size = ocp_nlp_reg_config_calculate_size();
	void * config_mem = malloc(config_size);
	ocp_nlp_reg_config *config = ocp_nlp_reg_config_assign(config_mem);

//	ocp_nlp_reg_noreg_config_initialize_default(config);
//	ocp_nlp_reg_mirror_config_initialize_default(config);
	ocp_nlp_reg_project_config_initialize_default(config);
//	ocp_nlp_reg_convexify_config_initialize_default(config);

    /************************************************
     * dims
     ************************************************/

	int N = 5;
	int nx_ = 4;
	int nu_ = 2;

	int *nx = malloc((N+1)*sizeof(int));
	int *nu = malloc((N+1)*sizeof(int));

	for(ii=0; ii<=N; ii++)
		nx[ii] = nx_;

	for(ii=0; ii<N; ii++)
		nu[ii] = nu_;
	nu[N] = 0;

	acados_size_t dims_size = config->dims_calculate_size(N);
	void * dims_mem = malloc(dims_size);
	ocp_nlp_reg_dims *dims = config->dims_assign(N, dims_mem);

	for(ii=0; ii<=N; ii++)
	{
		config->dims_set(config, dims, ii, "nx", nx+ii);
		config->dims_set(config, dims, ii, "nu", nu+ii);
	}

    /************************************************
     * opts
     ************************************************/

	acados_size_t opts_size = config->opts_calculate_size();
	void * opts_mem = malloc(opts_size);
	void *opts = config->opts_assign(opts_mem);

	config->opts_initialize_default(config, dims, opts);

	double epsilon = 1e-4;
	config->opts_set(config, opts, "epsilon", &epsilon);

    /************************************************
     * memory
     ************************************************/

	acados_size_t memory_size = config->memory_calculate_size(config, dims, opts);
	void * memory_mem = malloc(memory_size);
	void *memory = config->memory_assign(config, dims, opts, memory_mem);

	struct blasfeo_dmat *RSQrq = malloc((N+1)*sizeof(struct blasfeo_dmat));
	for(ii=0; ii<=N; ii++)
	{
		blasfeo_allocate_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], RSQrq+ii);
		for(jj=0; jj<nu[ii]; jj++)
			BLASFEO_DMATEL(RSQrq+ii, jj, jj) = 2*ii;
		for(jj=0; jj<nx[ii]; jj++)
			BLASFEO_DMATEL(RSQrq+ii, nu[ii]+jj, nu[ii]+jj) = -ii;
	}

//	config->memory_set(config, dims, memory, "RSQrq_ptr", RSQrq);
//	config->memory_set(config, dims, memory, "BAbt_ptr", NULL);
	config->memory_set_RSQrq_ptr(dims, RSQrq, memory);
	config->memory_set_BAbt_ptr(dims, NULL, memory);

    /************************************************
     * regularize function
     ************************************************/

	printf("\nbefore regularization\n\n");
	for(ii=0; ii<=N; ii++)
		blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], RSQrq+ii, 0, 0);

	// regularization
	config->regularize(config, dims, opts, memory);

	printf("\nafter regularization\n\n");
	for(ii=0; ii<=N; ii++)
		blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], RSQrq+ii, 0, 0);

    /************************************************
     * free memory & return
     ************************************************/

	free(nx);
	free(nu);
	free(config_mem);
	free(dims_mem);
	free(opts_mem);
	free(memory_mem);
	for(ii=0; ii<=N; ii++)
		blasfeo_free_dmat(RSQrq+ii);
	free(RSQrq);

	printf("\nsuccess !\n\n");

	return 0;

}
