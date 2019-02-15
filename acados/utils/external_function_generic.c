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

#include <assert.h>
#include <stdlib.h>

#include "acados/utils/external_function_generic.h"
#include "acados/utils/mem.h"

/************************************************
 * generic external function
 ************************************************/



/************************************************
 * casadi utils
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
            const int *idxcol = sparsity + 2;
            for (int i = 0; i < ncol; ++i)
            {
                nnz += idxcol[i + 1] - idxcol[i];
            }
        }
    }

    return nnz;
}



static void d_cvt_casadi_to_colmaj(double *in, int *sparsity_in, double *out)
{
    int ii, jj, idx;

    int nrow = sparsity_in[0];
    int ncol = sparsity_in[1];
    int dense = sparsity_in[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    if (dense)
    {
        for (ii = 0; ii < ncol * nrow; ii++) out[ii] = in[ii];
    }
    else
    {
        double *ptr = in;
        int *idxcol = sparsity_in + 2;
        int *row = sparsity_in + ncol + 3;
        // Fill with zeros
        for (jj = 0; jj < ncol; jj++)
            for (ii = 0; ii < nrow; ii++) out[ii + jj * nrow] = 0.0;
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                out[row[idx] + jj * nrow] = ptr[0];
                ptr++;
            }
        }
    }

    return;
}



static void d_cvt_colmaj_to_casadi(double *in, double *out, int *sparsity_out)
{
    int ii, jj, idx;

    int nrow = sparsity_out[0];
    int ncol = sparsity_out[1];
    int dense = sparsity_out[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    if (dense)
    {
        for (ii = 0; ii < ncol * nrow; ii++) out[ii] = in[ii];
    }
    else
    {
        double *ptr = out;
        int *idxcol = sparsity_out + 2;
        int *row = sparsity_out + ncol + 3;
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                ptr[0] = in[row[idx] + jj * nrow];
                ptr++;
            }
        }
    }

    return;
}



// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_casadi_to_dmat(double *in, int *sparsity_in, struct blasfeo_dmat *out)
{
    int jj, idx;

    int nrow = sparsity_in[0];
    int ncol = sparsity_in[1];
    int dense = sparsity_in[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    if (dense)
    {
        blasfeo_pack_dmat(nrow, ncol, in, nrow, out, 0, 0);
    }
    else
    {
        double *ptr = in;
        int *idxcol = sparsity_in + 2;
        int *row = sparsity_in + ncol + 3;
        // Fill with zeros
        blasfeo_dgese(nrow, ncol, 0.0, out, 0, 0);
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                BLASFEO_DMATEL(out, row[idx], jj) = ptr[0];
                ptr++;
            }
        }
    }

    return;
}



// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_dmat_to_casadi(struct blasfeo_dmat *in, double *out, int *sparsity_out)
{
    int jj, idx;

    int nrow = sparsity_out[0];
    int ncol = sparsity_out[1];
    int dense = sparsity_out[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    if (dense)
    {
        blasfeo_unpack_dmat(nrow, ncol, in, 0, 0, out, nrow);
    }
    else
    {
        double *ptr = out;
        int *idxcol = sparsity_out + 2;
        int *row = sparsity_out + ncol + 3;
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                ptr[0] = BLASFEO_DMATEL(in, row[idx], jj);
                ptr++;
            }
        }
    }

    return;
}



// column vector: assume sparsity_in[1] = 1 !!! or empty vector;
// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_casadi_to_dvec(double *in, int *sparsity_in, struct blasfeo_dvec *out)
{
    int idx;

    assert((sparsity_in[1] == 1) | (sparsity_in[0] == 0) | (sparsity_in[1] == 0));

    int n = sparsity_in[0];
    int dense = sparsity_in[2];

    if (n<=0)
        return;

    if (dense)
    {
        blasfeo_pack_dvec(n, in, out, 0);
    }
    else
    {
        double *ptr = in;
        int *idxcol = sparsity_in + 2;
        int *row = sparsity_in + 1 + 3;
        // Fill with zeros
        blasfeo_dvecse(n, 0.0, out, 0);
        // Copy nonzeros
        for (idx = idxcol[0]; idx != idxcol[1]; idx++)
        {
            BLASFEO_DVECEL(out, row[idx]) = ptr[0];
            ptr++;
        }
    }

    return;
}



// column vector: assume sparsity_in[1] = 1 !!! or empty vector;
// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_dvec_to_casadi(struct blasfeo_dvec *in, double *out, int *sparsity_out)
{
    int idx;

    assert((sparsity_out[1] == 1) | (sparsity_out[0] == 0) | (sparsity_out[1] == 0));

    int n = sparsity_out[0];
    int dense = sparsity_out[2];

    if (n<=0)
        return;

    if (dense)
    {
        blasfeo_unpack_dvec(n, in, 0, out);
    }
    else
    {
        double *ptr = out;
        int *idxcol = sparsity_out + 2;
        int *row = sparsity_out + 1 + 3;
        // Copy nonzeros
        for (idx = idxcol[0]; idx != idxcol[1]; idx++)
        {
            ptr[0] = BLASFEO_DVECEL(in, row[idx]);
            ptr++;
        }
    }

    return;
}



static void d_cvt_casadi_to_colmaj_args(double *in, int *sparsity_in, struct colmaj_args *out)
{
    int ii, jj, idx;

    int nrow = sparsity_in[0];
    int ncol = sparsity_in[1];
    int dense = sparsity_in[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    double *A = out->A;
    int lda = out->lda;

    if (dense)
    {
        for (ii = 0; ii < ncol; ii++)
            for (jj = 0; jj < nrow; jj++) A[ii + jj * lda] = in[ii + ncol * jj];
    }
    else
    {
        double *ptr = in;
        int *idxcol = sparsity_in + 2;
        int *row = sparsity_in + ncol + 3;
        // Fill with zeros
        for (jj = 0; jj < ncol; jj++)
            for (ii = 0; ii < nrow; ii++) A[ii + jj * lda] = 0.0;
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                A[row[idx] + jj * lda] = ptr[0];
                ptr++;
            }
        }
    }

    return;
}



static void d_cvt_colmaj_args_to_casadi(struct colmaj_args *in, double *out, int *sparsity_out)
{
    int ii, jj, idx;

    int nrow = sparsity_out[0];
    int ncol = sparsity_out[1];
    int dense = sparsity_out[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    double *A = in->A;
    int lda = in->lda;

    if (dense)
    {
        for (ii = 0; ii < ncol; ii++)
            for (jj = 0; jj < nrow; jj++) out[ii + ncol * jj] = A[ii + jj * lda];
    }
    else
    {
        double *ptr = out;
        int *idxcol = sparsity_out + 2;
        int *row = sparsity_out + ncol + 3;
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                ptr[0] = A[row[idx] + jj * lda];
                ptr++;
            }
        }
    }

    return;
}



// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_casadi_to_dmat_args(double *in, int *sparsity_in, struct blasfeo_dmat_args *out)
{
    int jj, idx;

    int nrow = sparsity_in[0];
    int ncol = sparsity_in[1];
    int dense = sparsity_in[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    struct blasfeo_dmat *A = out->A;
    int ai = out->ai;
    int aj = out->aj;

    if (dense)
    {
        blasfeo_pack_dmat(nrow, ncol, in, nrow, A, ai, aj);
    }
    else
    {
        double *ptr = in;
        int *idxcol = sparsity_in + 2;
        int *row = sparsity_in + ncol + 3;
        // Fill with zeros
        blasfeo_dgese(nrow, ncol, 0.0, A, ai, aj);
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                BLASFEO_DMATEL(A, ai + row[idx], aj + jj) = ptr[0];
                ptr++;
            }
        }
    }

    return;
}



// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_dmat_args_to_casadi(struct blasfeo_dmat_args *in, double *out, int *sparsity_out)
{
    int jj, idx;

    int nrow = sparsity_out[0];
    int ncol = sparsity_out[1];
    int dense = sparsity_out[2];

    if ((nrow<=0 )| (ncol<=0))
        return;

    struct blasfeo_dmat *A = in->A;
    int ai = in->ai;
    int aj = in->aj;

    if (dense)
    {
        blasfeo_unpack_dmat(nrow, ncol, A, ai, aj, out, nrow);
    }
    else
    {
        double *ptr = out;
        int *idxcol = sparsity_out + 2;
        int *row = sparsity_out + ncol + 3;
        // Copy nonzeros
        for (jj = 0; jj < ncol; jj++)
        {
            for (idx = idxcol[jj]; idx != idxcol[jj + 1]; idx++)
            {
                ptr[0] = BLASFEO_DMATEL(A, ai + row[idx], aj + nrow);
                ptr++;
            }
        }
    }

    return;
}



// column vector: assume sparsity_in[1] = 1 !!! or empty vector;
// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_casadi_to_dvec_args(double *in, int *sparsity_in, struct blasfeo_dvec_args *out)
{
    int idx;

    assert((sparsity_in[1] == 1) | (sparsity_in[0] == 0) | (sparsity_in[1] == 0));

    int n = sparsity_in[0];
    int dense = sparsity_in[2];

    if (n<=0)
        return;

    struct blasfeo_dvec *x = out->x;
    int xi = out->xi;

    if (dense)
    {
        blasfeo_pack_dvec(n, in, x, xi);
    }
    else
    {
        double *ptr = in;
        int *idxcol = sparsity_in + 2;
        int *row = sparsity_in + 1 + 3;
        // Fill with zeros
        blasfeo_dvecse(n, 0.0, x, xi);
        // Copy nonzeros
        for (idx = idxcol[0]; idx != idxcol[1]; idx++)
        {
            BLASFEO_DVECEL(x, xi + row[idx]) = ptr[0];
            ptr++;
        }
    }

    return;
}



// column vector: assume sparsity_in[1] = 1 !!! or empty vector;
// TODO(all): detect if dense from number of elements per column !!!
static void d_cvt_dvec_args_to_casadi(struct blasfeo_dvec_args *in, double *out, int *sparsity_out)
{
    int idx;

    assert((sparsity_out[1] == 1) | (sparsity_out[0] == 0) | (sparsity_out[1] == 0));

    int n = sparsity_out[0];
    int dense = sparsity_out[2];

    if (n<=0)
        return;

    struct blasfeo_dvec *x = in->x;
    int xi = in->xi;

    if (dense)
    {
        blasfeo_unpack_dvec(n, x, xi, out);
    }
    else
    {
        double *ptr = out;
        int *idxcol = sparsity_out + 2;
        int *row = sparsity_out + 1 + 3;
        // Copy nonzeros
        for (idx = idxcol[0]; idx != idxcol[1]; idx++)
        {
            ptr[0] = BLASFEO_DVECEL(x, xi + row[idx]);
            ptr++;
        }
    }

    return;
}



/************************************************
 * casadi external function
 ************************************************/

int external_function_casadi_struct_size()
{
    return sizeof(external_function_casadi);
}



void external_function_casadi_set_fun(external_function_casadi *fun, void *value)
{
    fun->casadi_fun = value;
    return;
}



void external_function_casadi_set_work(external_function_casadi *fun, void *value)
{
    fun->casadi_work = value;
    return;
}



void external_function_casadi_set_sparsity_in(external_function_casadi *fun, void *value)
{
    fun->casadi_sparsity_in = value;
    return;
}



void external_function_casadi_set_sparsity_out(external_function_casadi *fun, void *value)
{
    fun->casadi_sparsity_out = value;
    return;
}



void external_function_casadi_set_n_in(external_function_casadi *fun, void *value)
{
    fun->casadi_n_in = value;
    return;
}



void external_function_casadi_set_n_out(external_function_casadi *fun, void *value)
{
    fun->casadi_n_out = value;
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
    for (ii = 0; ii < fun->args_num; ii++)
        fun->args_size_tot += casadi_nnz(fun->casadi_sparsity_in(ii));

    // res
    fun->res_size_tot = 0;
    for (ii = 0; ii < fun->res_num; ii++)
        fun->res_size_tot += casadi_nnz(fun->casadi_sparsity_out(ii));

    int size = 0;

    // double pointers
    size += fun->args_num * sizeof(double *);  // args
    size += fun->res_num * sizeof(double *);   // res

    // ints
    size += fun->args_num * sizeof(int);  // args_size
    size += fun->res_num * sizeof(int);   // res_size
    size += fun->iw_size * sizeof(int);   // iw

    // doubles
    size += fun->args_size_tot * sizeof(double);  // args
    size += fun->res_size_tot * sizeof(double);   // res
    size += fun->w_size * sizeof(double);         // w

    size += 8;  // initial align
    size += 8;  // align to double

    // make_int_multiple_of(64, &size);

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
    assign_and_advance_double_ptrs(fun->args_num, &fun->args, &c_ptr);
    // res
    assign_and_advance_double_ptrs(fun->res_num, &fun->res, &c_ptr);

    // args_size
    assign_and_advance_int(fun->args_num, &fun->args_size, &c_ptr);
    for (ii = 0; ii < fun->args_num; ii++)
        fun->args_size[ii] = casadi_nnz(fun->casadi_sparsity_in(ii));
    // res_size
    assign_and_advance_int(fun->res_num, &fun->res_size, &c_ptr);
    for (ii = 0; ii < fun->res_num; ii++)
        fun->res_size[ii] = casadi_nnz(fun->casadi_sparsity_out(ii));
    // iw
    assign_and_advance_int(fun->iw_size, &fun->iw, &c_ptr);

    // align to double
    align_char_to(8, &c_ptr);

    // args
    for (ii = 0; ii < fun->args_num; ii++)
        assign_and_advance_double(fun->args_size[ii], &fun->args[ii], &c_ptr);
    // res
    for (ii = 0; ii < fun->res_num; ii++)
        assign_and_advance_double(fun->res_size[ii], &fun->res[ii], &c_ptr);
    // w
    assign_and_advance_double(fun->w_size, &fun->w, &c_ptr);

    assert((char *) raw_memory + external_function_casadi_calculate_size(fun) >= c_ptr);

    return;
}



void external_function_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                      ext_fun_arg_t *type_out, void **out)
{
    // cast into external casadi function
    external_function_casadi *fun = self;

    // loop index
    int ii;

    // in as args
    for (ii = 0; ii < fun->in_num; ii++)
    {
        switch (type_in[ii])
        {
            case COLMAJ:
                d_cvt_colmaj_to_casadi(in[ii], (double *) fun->args[ii],
                                       (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DMAT:
                d_cvt_dmat_to_casadi(in[ii], (double *) fun->args[ii],
                                     (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DVEC:
                d_cvt_dvec_to_casadi(in[ii], (double *) fun->args[ii],
                                     (int *) fun->casadi_sparsity_in(ii));
                break;
            case COLMAJ_ARGS:
                d_cvt_colmaj_args_to_casadi(in[ii], (double *) fun->args[ii],
                                            (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DMAT_ARGS:
                d_cvt_dmat_args_to_casadi(in[ii], (double *) fun->args[ii],
                                          (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DVEC_ARGS:
                d_cvt_dvec_args_to_casadi(in[ii], (double *) fun->args[ii],
                                          (int *) fun->casadi_sparsity_in(ii));
                break;

            case IGNORE_ARGUMENT:
                // do nothing
                break;

            default:
                printf("\ntype in %d\n", type_in[ii]);
                printf("\nUnknown external function argument type\n\n");
                exit(1);
        }
    }

    // call casadi function
    fun->casadi_fun((const double **) fun->args, fun->res, fun->iw, fun->w, NULL);

    for (ii = 0; ii < fun->out_num; ii++)
    {
        switch (type_out[ii])
        {
            case COLMAJ:
                d_cvt_casadi_to_colmaj((double *) fun->res[ii],
                                       (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case BLASFEO_DMAT:
                d_cvt_casadi_to_dmat((double *) fun->res[ii], (int *) fun->casadi_sparsity_out(ii),
                                     out[ii]);
                break;

            case BLASFEO_DVEC:
                d_cvt_casadi_to_dvec((double *) fun->res[ii], (int *) fun->casadi_sparsity_out(ii),
                                     out[ii]);
                break;
            case COLMAJ_ARGS:
                d_cvt_casadi_to_colmaj_args((double *) fun->res[ii],
                                            (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case BLASFEO_DMAT_ARGS:
                d_cvt_casadi_to_dmat_args((double *) fun->res[ii],
                                          (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case BLASFEO_DVEC_ARGS:
                d_cvt_casadi_to_dvec_args((double *) fun->res[ii],
                                          (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case IGNORE_ARGUMENT:
                // do nothing
                break;

            default:
                printf("\ntype out %d\n", type_out[ii]);
                printf("\nUnknown external function argument type\n\n");
                exit(1);
        }
    }

    return;
}

/************************************************
 * casadi external parametric function
 ************************************************/

int external_function_param_casadi_struct_size()
{
    return sizeof(external_function_param_casadi);
}



void external_function_param_casadi_set_fun(external_function_param_casadi *fun, void *value)
{
    fun->casadi_fun = value;
    return;
}



void external_function_param_casadi_set_work(external_function_param_casadi *fun, void *value)
{
    fun->casadi_work = value;
    return;
}



void external_function_param_casadi_set_sparsity_in(external_function_param_casadi *fun, void *value)
{
    fun->casadi_sparsity_in = value;
    return;
}



void external_function_param_casadi_set_sparsity_out(external_function_param_casadi *fun, void *value)
{
    fun->casadi_sparsity_out = value;
    return;
}



void external_function_param_casadi_set_n_in(external_function_param_casadi *fun, void *value)
{
    fun->casadi_n_in = value;
    return;
}



void external_function_param_casadi_set_n_out(external_function_param_casadi *fun, void *value)
{
    fun->casadi_n_out = value;
    return;
}



int external_function_param_casadi_calculate_size(external_function_param_casadi *fun, int np)
{
    // loop index
    int ii;

    // casadi wrapper as evaluate function
    fun->evaluate = &external_function_param_casadi_wrapper;

    // set param function
    fun->set_param = &external_function_param_casadi_set_param;

    // set number of parameters
    fun->np = np;

    fun->casadi_work(&fun->args_num, &fun->res_num, &fun->iw_size, &fun->w_size);

    fun->in_num = fun->casadi_n_in();
    fun->out_num = fun->casadi_n_out();

    // args
    fun->args_size_tot = 0;
    for (ii = 0; ii < fun->args_num; ii++)
        fun->args_size_tot += casadi_nnz(fun->casadi_sparsity_in(ii));

    // res
    fun->res_size_tot = 0;
    for (ii = 0; ii < fun->res_num; ii++)
        fun->res_size_tot += casadi_nnz(fun->casadi_sparsity_out(ii));

    int size = 0;

    // double pointers
    size += fun->args_num * sizeof(double *);  // args
    size += fun->res_num * sizeof(double *);   // res

    // ints
    size += fun->args_num * sizeof(int);  // args_size
    size += fun->res_num * sizeof(int);   // res_size
    size += fun->iw_size * sizeof(int);   // iw

    // doubles
    size += fun->args_size_tot * sizeof(double);  // args
    size += fun->res_size_tot * sizeof(double);   // res
    size += fun->w_size * sizeof(double);         // w
    size += fun->np * sizeof(double);             // p

    size += 8;  // initial align
    size += 8;  // align to double

    // make_int_multiple_of(64, &size);

    return size;
}



void external_function_param_casadi_assign(external_function_param_casadi *fun, void *raw_memory)
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
    assign_and_advance_double_ptrs(fun->args_num, &fun->args, &c_ptr);
    // res
    assign_and_advance_double_ptrs(fun->res_num, &fun->res, &c_ptr);

    // args_size
    assign_and_advance_int(fun->args_num, &fun->args_size, &c_ptr);
    for (ii = 0; ii < fun->args_num; ii++)
        fun->args_size[ii] = casadi_nnz(fun->casadi_sparsity_in(ii));
    // res_size
    assign_and_advance_int(fun->res_num, &fun->res_size, &c_ptr);
    for (ii = 0; ii < fun->res_num; ii++)
        fun->res_size[ii] = casadi_nnz(fun->casadi_sparsity_out(ii));
    // iw
    assign_and_advance_int(fun->iw_size, &fun->iw, &c_ptr);

    // align to double
    align_char_to(8, &c_ptr);

    // args
    for (ii = 0; ii < fun->args_num; ii++)
        assign_and_advance_double(fun->args_size[ii], &fun->args[ii], &c_ptr);
    // res
    for (ii = 0; ii < fun->res_num; ii++)
        assign_and_advance_double(fun->res_size[ii], &fun->res[ii], &c_ptr);
    // w
    assign_and_advance_double(fun->w_size, &fun->w, &c_ptr);
    // p
    assign_and_advance_double(fun->np, &fun->p, &c_ptr);

    assert((char *) raw_memory + external_function_param_casadi_calculate_size(fun, fun->np) >=
           c_ptr);

    return;
}



void external_function_param_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                            ext_fun_arg_t *type_out, void **out)
{
    // cast into external casadi function
    external_function_param_casadi *fun = self;

    // loop index
    int ii, jj;

    // in as args
    // skip last argument (that is the parameters vector)
    for (ii = 0; ii < fun->in_num - 1; ii++)
    {
        switch (type_in[ii])
        {
            case COLMAJ:
                d_cvt_colmaj_to_casadi(in[ii], (double *) fun->args[ii],
                                       (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DMAT:
                d_cvt_dmat_to_casadi(in[ii], (double *) fun->args[ii],
                                     (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DVEC:
                d_cvt_dvec_to_casadi(in[ii], (double *) fun->args[ii],
                                     (int *) fun->casadi_sparsity_in(ii));
                break;
            case COLMAJ_ARGS:
                d_cvt_colmaj_args_to_casadi(in[ii], (double *) fun->args[ii],
                                            (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DMAT_ARGS:
                d_cvt_dmat_args_to_casadi(in[ii], (double *) fun->args[ii],
                                          (int *) fun->casadi_sparsity_in(ii));
                break;

            case BLASFEO_DVEC_ARGS:
                d_cvt_dvec_args_to_casadi(in[ii], (double *) fun->args[ii],
                                          (int *) fun->casadi_sparsity_in(ii));
                break;

            default:
                printf("\nUnknown external function argument type\n\n");
                exit(1);
        }
    }
    // copy parametrs vector as last arg
    ii = fun->in_num - 1;
    for (jj = 0; jj < fun->np; jj++) fun->args[ii][jj] = fun->p[jj];

    // call casadi function
    fun->casadi_fun((const double **) fun->args, fun->res, fun->iw, fun->w, NULL);

    for (ii = 0; ii < fun->out_num; ii++)
    {
        switch (type_out[ii])
        {
            case COLMAJ:
                d_cvt_casadi_to_colmaj((double *) fun->res[ii],
                                       (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case BLASFEO_DMAT:
                d_cvt_casadi_to_dmat((double *) fun->res[ii], (int *) fun->casadi_sparsity_out(ii),
                                     out[ii]);
                break;

            case BLASFEO_DVEC:
                d_cvt_casadi_to_dvec((double *) fun->res[ii], (int *) fun->casadi_sparsity_out(ii),
                                     out[ii]);
                break;
            case COLMAJ_ARGS:
                d_cvt_casadi_to_colmaj_args((double *) fun->res[ii],
                                            (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case BLASFEO_DMAT_ARGS:
                d_cvt_casadi_to_dmat_args((double *) fun->res[ii],
                                          (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            case BLASFEO_DVEC_ARGS:
                d_cvt_casadi_to_dvec_args((double *) fun->res[ii],
                                          (int *) fun->casadi_sparsity_out(ii), out[ii]);
                break;

            default:
                printf("\nUnknown external function argument type\n\n");
                exit(1);
        }
    }

    return;
}



void external_function_param_casadi_set_param(void *self, double *p)
{
    // cast into external casadi function
    external_function_param_casadi *fun = self;

    // loop index
    int ii;

    for (ii = 0; ii < fun->np; ii++) fun->p[ii] = p[ii];

    return;
}
