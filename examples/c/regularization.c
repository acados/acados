// external
#include <stdio.h>
#include <stdlib.h>

// acados
#include <acados/ocp_nlp/ocp_nlp_reg_common.h>
#include <acados/ocp_nlp/ocp_nlp_reg_mirror.h>
#include <acados/ocp_nlp/ocp_nlp_reg_conv.h>



int main()
{
	
	printf("\nregularization example\n\n");

	int ii;

    /************************************************
     * config
     ************************************************/

	int config_size = ocp_nlp_reg_config_calculate_size();
	void * config_mem = malloc(config_size);
	ocp_nlp_reg_config *config = ocp_nlp_reg_config_assign(config_mem);

	ocp_nlp_reg_mirror_config_initialize_default(config);
//	ocp_nlp_reg_conv_config_initialize_default(config);

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

	int dims_size = config->dims_calculate_size(N);
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

	int opts_size = config->opts_calculate_size();
	void * opts_mem = malloc(opts_size);
	ocp_nlp_reg_opts *opts = config->opts_assign(opts_mem);

	config->opts_initialize_default(config, dims, opts);

	double delta = 1e-4;
	config->opts_set(config, dims, opts, "delta", &delta);



    /************************************************
     * free memory & return
     ************************************************/

	free(nx);
	free(nu);
	free(config_mem);
	free(dims_mem);

	printf("\nsuccess !\n\n");

	return 0;

}
