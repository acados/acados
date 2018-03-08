/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/utils/mem.h"
#include "acados/utils/external_function_generic.h"

#include <assert.h>

#include <stdlib.h>



/************************************************
* generic external function
************************************************/



/************************************************
* casadi external function
************************************************/

static int casadi_nnz(const int *sparsity)
{
    int nnz = 0;
    if (sparsity != NULL)
	{
        const int nrow = sparsity[0];
        const int ncol = sparsity[1];
        const int dense = sparsity[2];
        if (dense)
		{
            nnz = nrow * ncol;
        }
		else
		{
            const int *colind = sparsity + 2;
            for (int i = 0; i < ncol; ++i)
			{
                nnz += colind[i + 1] - colind[i];
            }
        }
    }

    return nnz;
}



static void casadi_densify(const double *sparse_in, double *dense_out, const int *sparsity)
{
    const int_t nrow = sparsity[0];
    const int_t ncol = sparsity[1];
    const int_t dense = sparsity[2];

    if (dense)
	{
        for (int_t i = 0; i < ncol * nrow; i++)
			dense_out[i] = sparse_in[i];
    }
	else
	{
        // Fill with zeros
        for (int i = 0; i < ncol; i++)
            for (int j = 0; j < nrow; j++)
				dense_out[i * nrow + j] = 0.0;
        // Additional data
        const double *x = sparse_in;
        const int *colind = sparsity + 2;
		const int *row = sparsity + ncol + 3;
        // Copy nonzeros
        for (int i = 0; i < ncol; ++i)
            for (int el = colind[i]; el != colind[i + 1]; ++el)
                dense_out[row[el] + i * nrow] = *x++;
    }

	return;
}



int external_function_casadi_calculate_size(external_function_casadi *fun)
{

	// casadi wrapper as evaluate
	fun->evaluate = &external_function_casadi_wrapper;

	// loop index
	int ii;

	fun->casadi_work(&fun->args_num, &fun->res_num, &fun->iw_size, &fun->w_size);

	fun->in_num = fun->casadi_n_in();
	fun->out_num = fun->casadi_n_out();

	// args
	fun->args_size_tot = 0;
	for (ii=0; ii<fun->args_num; ii++)
		fun->args_size_tot += casadi_nnz(fun->casadi_sparsity_in(ii));

	// res
	fun->res_size_tot = 0;
	for (ii=0; ii<fun->res_num; ii++)
		fun->res_size_tot += casadi_nnz(fun->casadi_sparsity_out(ii));

	int size = 0;

	// double pointers
	size += fun->args_num*sizeof(double *); // args
	size += fun->res_num*sizeof(double *); // res

	// ints
	size += fun->args_num*sizeof(int); // args_size
	size += fun->res_num*sizeof(int); // res_size
	size += fun->iw_size*sizeof(int); // iw

	// doubles
	size += fun->args_size_tot*sizeof(double); // args
	size += fun->res_size_tot*sizeof(double); // res
	size += fun->w_size*sizeof(double); // w

    size += 8; // initial align
    size += 8; // align to double

//	make_int_multiple_of(64, &size);

	return size;

}



void external_function_casadi_assign(external_function_casadi *fun, void *raw_memory)
{

	// loop index
	int ii;

	// save initial pointer to external memory
	fun->ptr_ext_mem = raw_memory;

	// char pointer for byte advances
	char *c_ptr = raw_memory;

	// double pointers

	// initial align
    align_char_to(8, &c_ptr);

	// args
	assign_double_ptrs(fun->args_num, &fun->args, &c_ptr);
	// res
	assign_double_ptrs(fun->res_num, &fun->res, &c_ptr);

	// args_size
	assign_int(fun->args_num, &fun->args_size, &c_ptr);
	for (ii=0; ii<fun->args_num; ii++)
		fun->args_size[ii] = casadi_nnz(fun->casadi_sparsity_in(ii));
	// res_size
	assign_int(fun->res_num, &fun->res_size, &c_ptr);
	for (ii=0; ii<fun->res_num; ii++)
		fun->res_size[ii] = casadi_nnz(fun->casadi_sparsity_out(ii));
	// iw
	assign_int(fun->iw_size, &fun->iw, &c_ptr);

	// align to double
    align_char_to(8, &c_ptr);

	// args
	for (ii=0; ii<fun->args_num; ii++)
		assign_double(fun->args_size[ii], &fun->args[ii], &c_ptr);
	// res
	for (ii=0; ii<fun->res_num; ii++)
		assign_double(fun->res_size[ii], &fun->res[ii], &c_ptr);
	// w
	assign_double(fun->w_size, &fun->w, &c_ptr);

    assert((char *) raw_memory + external_function_casadi_calculate_size(fun) >= c_ptr);

	return;

}



void external_function_casadi_wrapper(void *self, double *in, double *out)
{

	// cast into external casadi function
	external_function_casadi *fun = self;

	// loop index
	int ii, jj;

	char *c_ptr;

	// in as args
	// TODO implement sparsify of input instead of copying them
	double *d_ptr = in;
	for (ii=0; ii<fun->in_num; ii++)
	{
		for (jj=0; jj<fun->args_size[ii]; jj++)
			fun->args[ii][jj] = d_ptr[jj];
		d_ptr += fun->args_size[ii];
	}

	// call casadi function
	fun->casadi_fun((const double **)fun->args, fun->res, fun->iw, fun->w, 0);

	const int *sparsity;
	int nrow, ncol;
	double *ptr_out = out;
	for (ii=0; ii<fun->out_num; ii++)
	{
		sparsity = fun->casadi_sparsity_out(ii);
		nrow = sparsity[0];
		ncol = sparsity[1];
		casadi_densify(fun->res[ii], ptr_out, sparsity); 
		ptr_out += nrow*ncol;
	}

	return;

}



