/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/math.h"
#include "acados/utils/types.h"

#if defined(__MABX2__)
double fmax(double a, double b)
{
    return a > b ? a : b;
}

int isnan(double x)
{
    return x != x;
}
#endif

void dgemm_nn_3l(int m, int n, int k, double *A, int lda, double *B, int ldb, double *C, int ldc)
{
    int ii, jj, kk;

    for (jj = 0; jj < n; jj++)
    {
        for (ii = 0; ii < m; ii++)
        {
            C[ii + ldc * jj] = 0;
            for (kk = 0; kk < k; kk++)
            {
                C[ii + ldc * jj] += A[ii + lda * kk] * B[kk + ldb * jj];
            }
        }
    }

    return;
}

// void dgemv_n_3l(int m, int n, double *A, int lda, double *x, double *y)
// {
//     int ii, jj;

//     for (jj = 0; jj < n; jj++)
//     {
//         for (ii = 0; ii < m; ii++)
//         {
//             y[ii] += A[ii + lda * jj] * x[jj];
//         }
//     }

//     return;
// }

// void dgemv_t_3l(int m, int n, double *A, int lda, double *x, double *y)
// {
//     int ii, jj;

//     for (ii = 0; ii < n; ii++)
//     {
//         for (jj = 0; jj < m; jj++)
//         {
//             y[ii] += A[jj + lda * ii] * x[jj];
//         }
//     }

//     return;
// }

// void dcopy_3l(int n, double *x, int incx, double *y, int incy)
// {
//     int ii;

//     if (incx == 1 && incy == 1)
//     {
//         for (ii = 0; ii < n; ii++)
//         {
//             y[ii] = x[ii];
//         }
//     }
//     else
//     {
//         for (ii = 0; ii < n; ii++)
//         {
//             y[ii * incy] = x[ii * incx];
//         }
//     }

//     return;
// }

void daxpy_3l(int n, double da, double *dx, double *dy)
{
    int i;
    for (i = 0; i < n; i++)
    {
        dy[i] += da * dx[i];
    }
}

void dscal_3l(int n, double da, double *dx)
{
    int i;
    for (i = 0; i < n; i++)
    {
        dx[i] *= da;
    }
}

/************************************************
 Routine that copies a matrix
************************************************/
void dmcopy(int row, int col, double *A, int lda, double *B, int ldb)
{
    int i, j;
    for (j = 0; j < col; j++)
    {
        for (i = 0; i < row; i++)
        {
            B[i + j * ldb] = A[i + j * lda];
        }
    }
}

int idamax_3l(int n, double *x)
{
    if (n <= 0) return 0;
    if (n == 1) return 0;

    double dabs;
    double dmax = (x[0] > 0 ? x[0] : -x[0]);
    int idmax = 0;
    int jj;
    for (jj = 1; jj < n; jj++)
    {
        dabs = (x[jj] > 0 ? x[jj] : -x[jj]);
        if (dabs > dmax)
        {
            dmax = dabs;
            idmax = jj;
        }
    }

    return idmax;
}

void dswap_3l(int n, double *x, int incx, double *y, int incy)
{
    if (n <= 0) return;

    double temp;
    int jj;
    for (jj = 0; jj < n; jj++)
    {
        temp = x[0];
        x[0] = y[0];
        y[0] = temp;
        x += incx;
        y += incy;
    }
}

void dger_3l(int m, int n, double alpha, double *x, int incx, double *y, int incy, double *A,
             int lda)
{
    if (m == 0 || n == 0 || alpha == 0.0) return;

    int i, j;
    double *px, *py, temp;

    py = y;
    for (j = 0; j < n; j++)
    {
        temp = alpha * py[0];
        px = x;
        for (i = 0; i < m; i++)
        {
            A[i + lda * j] += px[0] * temp;
            px += incx;
        }
        py += incy;
    }

    return;
}

void dgetf2_3l(int m, int n, double *A, int lda, int *ipiv, int *info)
{
    if (m <= 0 || n <= 0) return;

    int i, j, jp;

    double Ajj;

    int size_min = (m < n ? m : n);

    // find the pivot and test for singularity
    for (j = 0; j < size_min; j++)
    {
        jp = j + idamax_3l(m - j, &A[j + lda * j]);
        ipiv[j] = jp;
        if (A[jp + lda * j] != 0)
        {
            // apply the interchange to columns 0:n-1
            if (jp != j)
            {
                dswap_3l(n, &A[j], lda, &A[jp], lda);
            }
            // compute elements j+1:m-1 of j-th column
            if (j < m - 1)
            {
                Ajj = A[j + lda * j];
                if ((Ajj > 0 ? Ajj : -Ajj) >= 2.22e-16)
                {
                    dscal_3l(m - j - 1, 1.0 / Ajj, &A[j + 1 + lda * j]);
                }
                else
                {
                    for (i = j + 1; i < m; i++)
                    {
                        A[i + lda * j] /= Ajj;
                    }
                }
            }
        }
        else if (*info == 0)
        {
            *info = j + 1;
        }

        if (j < size_min)
        {
            // update trailing submatrix
            dger_3l(m - j - 1, n - j - 1, -1.0, &A[j + 1 + lda * j], 1, &A[j + lda * (j + 1)], lda,
                    &A[j + 1 + lda * (j + 1)], lda);
        }
    }

    return;
}

void dlaswp_3l(int n, double *A, int lda, int k1, int k2, int *ipiv)
{
    int i, j, k, ix, ix0, i1, i2, n32, ip;
    double temp;

    ix0 = k1;
    i1 = k1;
    i2 = k2;

    n32 = (n / 32) * 32;
    if (n32 != 0)
    {
        for (j = 0; j < n32; j += 32)
        {
            ix = ix0;
            for (i = i1; i < i2; i++)
            {
                ip = ipiv[ix];
                if (ip != i)
                {
                    for (k = j; k < j + 32; k++)
                    {
                        temp = A[i + lda * k];
                        A[i + lda * k] = A[ip + lda * k];
                        A[ip + lda * k] = temp;
                    }
                }
                ix++;
            }
        }
    }
    if (n32 != n)
    {
        ix = ix0;
        for (i = i1; i < i2; i++)
        {
            ip = ipiv[ix];
            if (ip != i)
            {
                for (k = n32; k < n; k++)
                {
                    temp = A[i + lda * k];
                    A[i + lda * k] = A[ip + lda * k];
                    A[ip + lda * k] = temp;
                }
            }
            ix++;
        }
    }

    return;
}

// left lower no-transp unit
void dtrsm_llnu_3l(int m, int n, double *A, int lda, double *B, int ldb)
{
    if (m == 0 || n == 0) return;

    int i, j, k;

    for (j = 0; j < n; j++)
    {
        for (k = 0; k < m; k++)
        {
            for (i = k + 1; i < m; i++)
            {
                B[i + ldb * j] -= B[k + ldb * j] * A[i + lda * k];
            }
        }
    }

    return;
}

// left upper no-transp non-unit
void dtrsm_lunn_3l(int m, int n, double *A, int lda, double *B, int ldb)
{
    if (m == 0 || n == 0) return;

    int i, j, k;

    for (j = 0; j < n; j++)
    {
        for (k = m - 1; k >= 0; k--)
        {
            B[k + ldb * j] /= A[k + lda * k];
            for (i = 0; i < k; i++)
            {
                B[i + ldb * j] -= B[k + ldb * j] * A[i + lda * k];
            }
        }
    }

    return;
}

void dgetrs_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb)
{
    if (n == 0 || nrhs == 0) return;

    // solve A * X = B

    // apply row interchanges to the rhs
    dlaswp_3l(nrhs, B, ldb, 0, n, ipiv);

    // solve L*X = B, overwriting B with X
    dtrsm_llnu_3l(n, nrhs, A, lda, B, ldb);

    // solve U*X = B, overwriting B with X
    dtrsm_lunn_3l(n, nrhs, A, lda, B, ldb);

    return;
}

void dgesv_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb, int *info)
{
    // compute the LU factorization of A
    dgetf2_3l(n, n, A, lda, ipiv, info);

    if (*info == 0)
    {
        // solve the system A*X = B, overwriting B with X
        dgetrs_3l(n, nrhs, A, lda, ipiv, B, ldb);
    }

    return;
}

/* one norm of a matrix */
double onenorm(int row, int col, double *ptrA)
{
    double max = 0.0;
    double temp;
    int i, j;
    temp = 0;
    for (j = 0; j < col; j++)
    {
        temp = fabs(*(ptrA + j * row));
        for (i = 1; i < row; i++)
        {
            temp += fabs(*(ptrA + j * row + i));
        }
        if (j == 0)
            max = temp;
        else if (temp > max)
            max = temp;
    }
    return max;
}

/* two norm of a vector */
// double twonormv(int n, double *ptrv)
// {
//     double temp;
//     temp = 0;
//     for (int i = 0; i < n; i++) temp += ptrv[i] * ptrv[i];
//     return (double) sqrt(temp);
// }

/* computes the Pade approximation of degree m of the matrix A */
void padeapprox(int m, int row, double *A)
{
    int row2 = row * row;
    /*    int i1 = 1;*/
    /*    double d0 = 0;*/
    /*    double d1 = 1;*/
    /*    double dm1 = -1;*/

    double *U = malloc(row * row * sizeof(double));
    for (int ii = 0; ii < row * row; ii++) U[ii] = 0.0;
    double *V = malloc(row * row * sizeof(double));
    for (int ii = 0; ii < row * row; ii++) V[ii] = 0.0;

    if (m == 3)
    {
        double c[] = {120, 60, 12, 1};
        double *A0 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A0[ii] = 0.0;
        for (int ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double *A2 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A2[ii] = 0.0;
        //        char ta = 'n'; double alpha = 1; double beta = 0;
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row,
        //        &beta, A2, &row);
        dgemm_nn_3l(row, row, row, A, row, A, row, A2, row);
        double *temp = malloc(row * row * sizeof(double));
        //        dscal_(&row2, &d0, temp, &i1);
        dscal_3l(row2, 0, temp);
        //        daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        daxpy_3l(row2, c[3], A2, temp);
        //        daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        daxpy_3l(row2, c[1], A0, temp);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp,
        //        &row, &beta, U, &row);
        dgemm_nn_3l(row, row, row, A, row, temp, row, U, row);
        //        dscal_(&row2, &d0, V, &i1);
        dscal_3l(row2, 0, V);
        //        daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        daxpy_3l(row2, c[2], A2, V);
        //        daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        daxpy_3l(row2, c[0], A0, V);
        free(A0);
        free(A2);
        free(temp);
    }
    else if (m == 5)
    {
        double c[] = {30240, 15120, 3360, 420, 30, 1};
        double *A0 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A0[ii] = 0.0;
        for (int ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double *A2 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A2[ii] = 0.0;
        double *A4 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A4[ii] = 0.0;
        //        char ta = 'n'; double alpha = 1; double beta = 0;
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row,
        //        &beta, A2, &row);
        dgemm_nn_3l(row, row, row, A, row, A, row, A2, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row,
        //        &beta, A4, &row);
        dgemm_nn_3l(row, row, row, A2, row, A2, row, A4, row);
        dmcopy(row, row, A4, row, V, row);
        double *temp = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) temp[ii] = 0.0;
        dmcopy(row, row, A4, row, temp, row);
        //        daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        daxpy_3l(row2, c[3], A2, temp);
        //        daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        daxpy_3l(row2, c[1], A0, temp);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp,
        //        &row, &beta, U, &row);
        dgemm_nn_3l(row, row, row, A, row, temp, row, U, row);
        //        dscal_(&row2, &c[4], V, &i1);
        dscal_3l(row2, c[4], V);
        //        daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        daxpy_3l(row2, c[2], A2, V);
        //        daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        daxpy_3l(row2, c[0], A0, V);
        free(A0);
        free(A2);
        free(A4);
        free(temp);
    }
    else if (m == 7)
    {
        double c[] = {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1};
        double *A0 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A0[ii] = 0.0;
        for (int ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double *A2 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A2[ii] = 0.0;
        double *A4 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A4[ii] = 0.0;
        double *A6 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A6[ii] = 0.0;
        //        char ta = 'n'; double alpha = 1; double beta = 1;
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row,
        //        &beta, A2, &row);
        dgemm_nn_3l(row, row, row, A, row, A, row, A2, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row,
        //        &beta, A4, &row);
        dgemm_nn_3l(row, row, row, A2, row, A2, row, A4, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A4, &row, A2, &row,
        //        &beta, A6, &row);
        dgemm_nn_3l(row, row, row, A4, row, A2, row, A6, row);
        double *temp = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) temp[ii] = 0.0;
        //        dscal_(&row2, &d0, temp, &i1);
        dscal_3l(row2, 0, temp);
        //        daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        daxpy_3l(row2, c[3], A2, temp);
        //        daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        daxpy_3l(row2, c[1], A0, temp);
        //        daxpy_(&row2, &c[5], A4, &i1, temp, &i1);
        daxpy_3l(row2, c[5], A4, temp);
        //        daxpy_(&row2, &c[7], A6, &i1, temp, &i1);
        daxpy_3l(row2, c[7], A6, temp);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp,
        //        &row, &beta, U, &row);
        dgemm_nn_3l(row, row, row, A, row, temp, row, U, row);
        //        dscal_(&row2, &d0, V, &i1);
        dscal_3l(row2, 0, V);
        //        daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        daxpy_3l(row2, c[2], A2, V);
        //        daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        daxpy_3l(row2, c[0], A0, V);
        //        daxpy_(&row2, &c[4], A4, &i1, V, &i1);
        daxpy_3l(row2, c[4], A4, V);
        //        daxpy_(&row2, &c[6], A6, &i1, V, &i1);
        daxpy_3l(row2, c[6], A6, V);
        free(A0);
        free(A2);
        free(A4);
        free(A6);
        free(temp);
    }
    else if (m == 9)
    {
        double c[] = {17643225600, 8821612800, 2075673600, 302702400, 30270240,
                      2162160,     110880,     3960,       90,        1};
        double *A0 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A0[ii] = 0.0;
        for (int ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double *A2 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A2[ii] = 0.0;
        double *A4 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A4[ii] = 0.0;
        double *A6 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A6[ii] = 0.0;
        double *A8 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A8[ii] = 0.0;
        //        char ta = 'n'; double alpha = 1; double beta = 0;
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row,
        //        &beta, A2, &row);
        dgemm_nn_3l(row, row, row, A, row, A, row, A2, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row,
        //        &beta, A4, &row);
        dgemm_nn_3l(row, row, row, A2, row, A2, row, A4, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A4, &row, A2, &row,
        //        &beta, A6, &row);
        dgemm_nn_3l(row, row, row, A4, row, A2, row, A6, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A6, &row, A2, &row,
        //        &beta, A8, &row);
        dgemm_nn_3l(row, row, row, A6, row, A2, row, A8, row);
        dmcopy(row, row, A8, row, V, row);
        double *temp = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) temp[ii] = 0.0;
        dmcopy(row, row, A8, row, temp, row);
        //        daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        daxpy_3l(row2, c[3], A2, temp);
        //        daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        daxpy_3l(row2, c[1], A0, temp);
        //        daxpy_(&row2, &c[5], A4, &i1, temp, &i1);
        daxpy_3l(row2, c[5], A4, temp);
        //        daxpy_(&row2, &c[7], A6, &i1, temp, &i1);
        daxpy_3l(row2, c[7], A6, temp);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp,
        //        &row, &beta, U, &row);
        dgemm_nn_3l(row, row, row, A, row, temp, row, U, row);
        //        dscal_(&row2, &c[8], V, &i1);
        dscal_3l(row2, c[8], V);
        //        daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        daxpy_3l(row2, c[2], A2, V);
        //        daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        daxpy_3l(row2, c[0], A0, V);
        //        daxpy_(&row2, &c[4], A4, &i1, V, &i1);
        daxpy_3l(row2, c[4], A4, V);
        //        daxpy_(&row2, &c[6], A6, &i1, V, &i1);
        daxpy_3l(row2, c[6], A6, V);
        free(A0);
        free(A2);
        free(A4);
        free(A6);
        free(A8);
        free(temp);
    }
    else if (m == 13)
    {  // tested
        double c[] = {64764752532480000,
                      32382376266240000,
                      7771770303897600,
                      1187353796428800,
                      129060195264000,
                      10559470521600,
                      670442572800,
                      33522128640,
                      1323241920,
                      40840800,
                      960960,
                      16380,
                      182,
                      1};
        double *A0 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A0[ii] = 0.0;
        for (int ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double *A2 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A2[ii] = 0.0;
        double *A4 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A4[ii] = 0.0;
        double *A6 = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) A6[ii] = 0.0;
        //        char ta = 'n'; double alpha = 1; double beta = 0;
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row,
        //        &beta, A2, &row);
        dgemm_nn_3l(row, row, row, A, row, A, row, A2, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row,
        //        &beta, A4, &row);
        dgemm_nn_3l(row, row, row, A2, row, A2, row, A4, row);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A4, &row, A2, &row,
        //        &beta, A6, &row);
        dgemm_nn_3l(row, row, row, A4, row, A2, row, A6, row);
        dmcopy(row, row, A2, row, U, row);
        double *temp = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) temp[ii] = 0.0;
        //        dscal_(&row2, &c[9], U, &i1);
        dscal_3l(row2, c[9], U);
        //        daxpy_(&row2, &c[11], A4, &i1, U, &i1);
        daxpy_3l(row2, c[11], A4, U);
        //        daxpy_(&row2, &c[13], A6, &i1, U, &i1);
        daxpy_3l(row2, c[13], A6, U);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A6, &row, U, &row,
        //        &beta, temp, &row);
        dgemm_nn_3l(row, row, row, A6, row, U, row, temp, row);
        //        daxpy_(&row2, &c[7], A6, &i1, temp, &i1);
        daxpy_3l(row2, c[7], A6, temp);
        //        daxpy_(&row2, &c[5], A4, &i1, temp, &i1);
        daxpy_3l(row2, c[5], A4, temp);
        //        daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        daxpy_3l(row2, c[3], A2, temp);
        //        daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        daxpy_3l(row2, c[1], A0, temp);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp,
        //        &row, &beta, U, &row);
        dgemm_nn_3l(row, row, row, A, row, temp, row, U, row);
        dmcopy(row, row, A2, row, temp, row);
        //        dscal_(&row2, &c[8], V, &i1);
        dscal_3l(row2, c[8], V);
        //        daxpy_(&row2, &c[12], A6, &i1, temp, &i1);
        daxpy_3l(row2, c[12], A6, temp);
        //        daxpy_(&row2, &c[10], A4, &i1, temp, &i1);
        daxpy_3l(row2, c[10], A4, temp);
        //        dgemm_(&ta, &ta, &row, &row, &row, &alpha, A6, &row, temp,
        //        &row, &beta, V, &row);
        dgemm_nn_3l(row, row, row, A6, row, temp, row, V, row);
        //        daxpy_(&row2, &c[6], A6, &i1, V, &i1);
        daxpy_3l(row2, c[6], A6, V);
        //        daxpy_(&row2, &c[4], A4, &i1, V, &i1);
        daxpy_3l(row2, c[4], A4, V);
        //        daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        daxpy_3l(row2, c[2], A2, V);
        //        daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        daxpy_3l(row2, c[0], A0, V);
        free(A0);
        free(A2);
        free(A4);
        free(A6);
        free(temp);
    }
    else
    {
        printf("%s\n", "Wrong Pade approximatin degree");
        exit(1);
    }
    double *D = malloc(row * row * sizeof(double));
    for (int ii = 0; ii < row * row; ii++) D[ii] = 0.0;
    //    dcopy_(&row2, V, &i1, A, &i1);
    dmcopy(row, row, V, row, A, row);
    //    daxpy_(&row2, &d1, U, &i1, A, &i1);
    daxpy_3l(row2, 1.0, U, A);
    //    dcopy_(&row2, V, &i1, D, &i1);
    dmcopy(row, row, V, row, D, row);
    //    daxpy_(&row2, &dm1, U, &i1, D, &i1);
    daxpy_3l(row2, -1.0, U, D);
    int *ipiv = (int *) calloc(row, sizeof(int));
    int info = 0;
    //    dgesv_(&row, &row, D, &row, ipiv, A, &row, &info);
    dgesv_3l(row, row, D, row, ipiv, A, row, &info);
    free(ipiv);
    free(D);
    free(U);
    free(V);
}

void expm(int row, double *A)
{
    int i;

    int m_vals[] = {3, 5, 7, 9, 13};
    double theta[] = {0.01495585217958292, 0.2539398330063230, 0.9504178996162932,
                      2.097847961257068, 5.371920351148152};
    int lentheta = 5;

    double normA = onenorm(row, row, A);

    if (normA <= theta[4])
    {
        for (i = 0; i < lentheta; i++)
        {
            if (normA <= theta[i])
            {
                padeapprox(m_vals[i], row, A);
                break;
            }
        }
    }
    else
    {
        int s;
        double t = frexp(normA / (theta[4]), &s);
        s = s - (t == 0.5);
        t = pow(2, -s);
        int row2 = row * row;
        /*        int i1 = 1;*/
        //        dscal_(&row2, &t, A, &i1);
        dscal_3l(row2, t, A);
        padeapprox(m_vals[4], row, A);
        double *temp = malloc(row * row * sizeof(double));
        for (int ii = 0; ii < row * row; ii++) temp[ii] = 0.0;
        //        char ta = 'n'; double alpha = 1; double beta = 0;
        for (i = 0; i < s; i++)
        {
            //            dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A,
            //            &row, &beta, temp, &row);
            dgemm_nn_3l(row, row, row, A, row, A, row, temp, row);
            dmcopy(row, row, temp, row, A, row);
        }
        free(temp);
    }
}

// TODO(dimitris): move this to condensing module once implemented
// NOTE(oj): probably done!
// compute the memory size of condensing for [x u] order of bounds (instead of [u x] in hpipm)
// void d_compute_qp_size_ocp2dense_rev(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng,
//                                      int *nvd, int *ned, int *nbd, int *ngd)
// {
//     int ii, jj;

//     *nvd = 0;
//     *ned = 0;
//     *nbd = 0;
//     *ngd = 0;

//     // first stage
//     *nvd += nx[0] + nu[0];
//     *nbd += nb[0];
//     *ngd += ng[0];
//     // remaining stages
//     for (ii = 1; ii <= N; ii++)
//     {
//         *nvd += nu[ii];
//         for (jj = 0; jj < nb[ii]; jj++)
//         {
//             if (hidxb[ii][jj] < nx[ii])
//             {  // state constraint
//                 (*ngd)++;
//             }
//             else
//             {  // input constraint
//                 (*nbd)++;
//             }
//         }
//         *ngd += ng[ii];
//     }
// }



static double hypot2(double x, double y)
{
    return sqrt(x*x + y*y);
}



/* Symmetric Householder reduction to tridiagonal form. */
static void tred2(int dim, double *V, double *d, double *e)
{
    /* This is derived from the Algol procedures tred2 by
    Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    Fortran subroutine in EISPACK. */

    int i, j, k;
    double f, g, h, hh, scale;
    for (j = 0; j < dim; j++)
    {
        d[j] = V[(dim - 1) * dim + j];
    }

    /* Householder reduction to tridiagonal form. */

    for (i = dim - 1; i > 0; i--)
    {
        /* Scale to avoid under/overflow. */

        scale = 0.0;
        h = 0.0;
        for (k = 0; k < i; k++)
        {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0)
        {
            e[i] = d[i - 1];
            for (j = 0; j < i; j++)
            {
                d[j] = V[(i - 1) * dim + j];
                V[i * dim + j] = 0.0;
                V[j * dim + i] = 0.0;
            }
        }
        else
        {
            /* Generate Householder vector. */

            for (k = 0; k < i; k++)
            {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            f = d[i - 1];
            g = sqrt(h);
            if (f > 0)
            {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i - 1] = f - g;
            for (j = 0; j < i; j++)
            {
                e[j] = 0.0;
            }

            /* Apply similarity transformation to remaining columns. */

            for (j = 0; j < i; j++)
            {
                f = d[j];
                V[j * dim + i] = f;
                g = e[j] + V[j * dim + j] * f;
                for (k = j + 1; k <= i - 1; k++)
                {
                    g += V[k * dim + j] * d[k];
                    e[k] += V[k * dim + j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++)
            {
                e[j] /= h;
                f += e[j] * d[j];
            }
            hh = f / (h + h);
            for (j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }
            for (j = 0; j < i; j++)
            {
                f = d[j];
                g = e[j];
                for (k = j; k <= i - 1; k++)
                {
                    V[k * dim + j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[(i - 1) * dim + j];
                V[i * dim + j] = 0.0;
            }
        }
        d[i] = h;
    }

    /* Accumulate transformations. */

    for (i = 0; i < dim - 1; i++)
    {
        V[(dim - 1) * dim + i] = V[i * dim + i];
        V[i * dim + i] = 1.0;
        h = d[i + 1];
        if (h != 0.0)
        {
            for (k = 0; k <= i; k++)
            {
                d[k] = V[k * dim + i + 1] / h;
            }
            for (j = 0; j <= i; j++)
            {
                g = 0.0;
                for (k = 0; k <= i; k++)
                {
                    g += V[k * dim + i + 1] * V[k * dim + j];
                }
                for (k = 0; k <= i; k++)
                {
                    V[k * dim + j] -= g * d[k];
                }
            }
        }
        for (k = 0; k <= i; k++)
        {
            V[k * dim + i + 1] = 0.0;
        }
    }
    for (j = 0; j < dim; j++)
    {
        d[j] = V[(dim - 1) * dim + j];
        V[(dim - 1) * dim + j] = 0.0;
    }
    if (dim > 0)
    {
        V[(dim - 1) * dim + dim - 1] = 1.0;
        e[0] = 0.0;
    }
}



/* Symmetric tridiagonal QL algorithm. */
static void tql2(int dim, double *V, double *d, double *e)
{
    /*  This is derived from the Algol procedures tql2, by
    Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    Fortran subroutine in EISPACK. */
    // http://www.netlib.org/eispack/tql2.f

    int i, m, l, k;
    double g, p, r, dl1, h, f, tst1, eps;
    double c, c2, c3, el1, s, s2;

    for (i = 1; i < dim; i++)
    {
        e[i - 1] = e[i];
    }
    e[dim - 1] = 0.0;

    f = 0.0;
    tst1 = 0.0;
    eps = pow(2.0, -52.0);
    for (l = 0; l < dim; l++)
    {
        /* Find small subdiagonal element */

        tst1 = fmax(tst1, fabs(d[l]) + fabs(e[l]));
        m = l;
        while (m < dim-1)
        {
            if (fabs(e[m]) <= eps * tst1)
            {
                break;
            }
            m++;
        }

        /* If m == l, d[l] is an eigenvalue,
        otherwise, iterate. */

        if (m > l)
        {
            int iter = 0;
            do
            {
                iter = iter + 1;
                /* Compute implicit shift */

                g = d[l];
                p = (d[l + 1] - g) / (2.0 * e[l]);
                r = hypot2(p, 1.0);
                if (p < 0)
                {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                dl1 = d[l + 1];
                h = g - d[l];
                for (i = l + 2; i < dim; i++)
                {
                    d[i] -= h;
                }
                f = f + h;

                /* Implicit QL transformation. */

                p = d[m];
                c = 1.0;
                c2 = c;
                c3 = c;
                el1 = e[l + 1];
                s = 0.0;
                s2 = 0.0;
                for (i = m - 1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    /* Accumulate transformation. */

                    for (k = 0; k < dim; k++)
                    {
                        h = V[k * dim + i + 1];
                        V[k * dim + i + 1] = s * V[k * dim + i] + c * h;
                        V[k * dim + i] = c * V[k * dim + i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                /* Check for convergence. */

            } while (fabs(e[l]) > eps * tst1 && iter < 20); /* (Check iteration count here.) */
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }
}



void acados_eigen_decomposition(int dim, double *A, double *V, double *d, double *e)
{
    int i, j;

    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++)
            V[i*dim+j] = A[i*dim+j];

    tred2(dim, V, d, e);
    tql2(dim, V, d, e);

    return;
}



double minimum_of_doubles(double *x, int n)
{
    double min = x[0];
    for (int c = 1; c < n; c++)
    {
        if (x[c] < min)
        {
            min = x[c];
        }
    }
    return min;
}

void neville_algorithm(double xx, int n, double *x, double *Q, double *out)
{  // Neville's algorithm
    // writes value of interpolating polynom corresponding to the nodes (x_i, Q_i), i = 0,...,n
    // evaluated evaluated at xx into out
    for (int i = n; i > 0; i--)
    {
        for (int j = 0; j < i; j++)
        {
            Q[j] = (xx - x[j]) * Q[j + 1] -
                   (xx - x[j + n - i + 1]) * Q[j];
            Q[j] = Q[j] / (x[j + n - i + 1] - x[j]);
        }
    }
    out[0] = Q[0];
}
