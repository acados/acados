/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <qpOASES_e.h>


USING_NAMESPACE_QPOASES

int sumOfSquares()
{
	int i;

	/* sum of first n squares */
	const int N = 100;
	real_t av[N];
	real_t aT[N];
	real_t bv[N];
	real_t c;

	for (i = 0; i < N; i++) av[i] = i+1;
	for (i = 0; i < N; i++) aTv[i] = i+1;
	for (i = 0; i < N; i++) bv[i] = i+1;

	myStatic DenseMatrix a;
	DenseMatrixCON( &a, 1, N, N, av);
	myStatic DenseMatrix aT;
	DenseMatrixCON( &aT, N, 1, 1, aTv);

	a.times(1, 1.0, bv, N, 0.0, &c, 1);
	fprintf(stdFile, "Dot product; Error in sum of first %d squares: %9.2e\n", N, c - (1.0/6.0)*N*(N+1)*(2*N+1));
	aT.transTimes(1, 1.0, bv, N, 0.0, &c, 1);
	fprintf(stdFile, "Transpose; Error in sum of first %d squares: %9.2e\n", N, c - (1.0/6.0)*N*(N+1)*(2*N+1));

	return 0;
}

int hilbert()
{
	int i, j;
	real_t d, err;

	/* permuted 4x4 Hilbert matrix, row major format */
	real_t _Av[] = {1.0/3.0, 1.0, 0.5, 0.25,
		0.25, 0.5, 1.0/3.0, 0.2,
		0.2, 1.0/3.0, 0.25, 1.0/6.0,
		1.0/6.0, 0.25, 0.2, 1.0/7.0};
	/* and its inverse, column major format */
	real_t Bv[] = {240, 16, -120, -140,
		-2700, -120, 1200, 1680,
		6480, 240, -2700, -4200,
		-4200, -140, 1680, 2800};

	/* result */
	real_t *Av = new real_t[4*4];
	real_t *Cv = new real_t[4*4];

	myStatic DenseMatrix A;
	DenseMatrixCON( &A, 4, 4, 4, Av );

	for (i = 0; i < 16; i++) Av[i] = _Av[i];

	A.times(4, 1.0, Bv, 4, 0.0, Cv, 4);

	err = 0.0;
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 4; i++)
		{
			d = fabs(Cv[j*4+i] - static_cast<real_t>(i == j));
			if (d > err) err = d;
		}
	}
	fprintf(stdFile, "Hilbert; Deviation from identity: %9.2e\n", err);

	delete[] Cv;
	A.free ();	// or delete[] Av;

	return 0;
}

int submatrix()
{
	int i, j;
	real_t d, err;

	/* 2x3 transposed submatrix */
	real_t _Asubv[] = {1.0/3.0, 0.25,
		1.0, 0.5,
		0.5, 1.0/3.0,
		0.25, 0.2};
	real_t Bsubv[] = {240, 16, -120, -140,
		-2700, -120, 1200, 1680};
	real_t Csubv[2*2];

	real_t *Asubv = new real_t[2*4];

	for (i = 0; i < 8; i++) Asubv[i] = _Asubv[i];

	myStatic DenseMatrix Asub;
	DenseMatrixCON( &Asub, 4, 2, 2, Asubv );

	Asub.transTimes(2, 1.0, Bsubv, 4, 0.0, Csubv, 2);

	err = 0.0;
	for (j = 0; j < 2; j++)
	{
		for (i = 0; i < 2; i++)
		{
			d = fabs(Csubv[j*2+i] - static_cast<real_t>(i == j));
			if (d > err) err = d;
		}
	}
	fprintf(stdFile, "Submatrix transpose; Deviation from identity: %9.2e\n", err);
	
	Asub.free ();	// or delete[] Asubv;

	return 0;
}

int indexDenseSubmatrix()
{
	/* dense submatrices via index lists */
	const int M = 20, N = 15, K = 5;
	int i, j, k, m, n;
	Indexlist rows(N), cols(N);
	int *rNum, *cNum;
	real_t err = 0.0;
	real_t *Av, *X, *Y, *x, *y, *xc, *yc;

	// prepare index lists
	m = (M-3)/4;
	n = (N-3)/4;
	for (i = 0; i < m; i++) { rows.addNumber(M-1 - 2*i); rows.addNumber(2*i); }
	for (i = 0; i < n; i++) { cols.addNumber(N-1 - 2*i); cols.addNumber(2*i); }
	m *= 2;
	n *= 2;

	rows.getNumberArray(&rNum);
	fprintf(stdFile, "Rows: ");
	for (i = 0; i < m; i++) fprintf(stdFile, " %2d", rNum[i]);
	fprintf(stdFile, "\n");

	cols.getNumberArray(&cNum);
	fprintf(stdFile, "Cols: ");
	for (i = 0; i < n; i++) fprintf(stdFile, " %2d", cNum[i]);
	fprintf(stdFile, "\n");

	// prepare input matrices
	Av = new real_t[M*N];
	X = new real_t[n*K];
	Y = new real_t[m*K];
	x = new real_t[n*K];
	y = new real_t[m*K];
	xc = new real_t[n*K];
	yc = new real_t[m*K];

	myStatic DenseMatrix A;
	DenseMatrixCON( &A, M, N, N, Av );
	for (i = 0; i < M*N; i++) Av[i] = -0.5*N*M + (real_t)i;
	for (i = 0; i < n*K; i++) X[i] = 1.0 / (real_t)(i+1);
	for (i = 0; i < m*K; i++) Y[i] = 1.0 / (real_t)(i+1);

	// multiply
	A.times(&rows, &cols, K, 1.0, X, n, 0.0, y, m);

	// check result
	for (j = 0; j < m; j++)
	{
		for (k = 0; k < K; k++)
		{
			yc[j+k*m] = -y[j+k*m];
			for (i = 0; i < n; i++)
				yc[j+k*m] += Av[cNum[i]+rNum[j]*N] * X[i+k*n];
			if (fabs(yc[j+k*m]) > err) err = fabs(yc[j+k*m]);
		}
	}
	fprintf(stdFile, "Indexlist submatrix; error: %9.2e\n", err);

	// transpose multiply
	A.transTimes(&rows, &cols, K, 1.0, Y, m, 0.0, x, n);

	// check result
	err = 0.0;
	for (j = 0; j < n; j++)
	{
		for (k = 0; k < K; k++)
		{
			xc[j+k*n] = -x[j+k*n];
			for (i = 0; i < m; i++)
				xc[j+k*n] += Av[cNum[j]+rNum[i]*N] * Y[i+k*m];
			if (fabs(xc[j+k*n]) > err) err = fabs(xc[j+k*n]);
		}
	}
	fprintf(stdFile, "Indexlist transpose submatrix; error: %9.2e\n", err);

	// check result

	// clean up
	delete[] yc;
	delete[] xc;
	delete[] y;
	delete[] x;
	delete[] Y;
	delete[] X;
	A.free ();	// or delete[] Av;

	return 0;
}

long spGetCol()
{
	long j, i;
	sparse_int_t *ir = new sparse_int_t[10];
	sparse_int_t *jc = new sparse_int_t[4];
	real_t *val = new real_t[10];
	real_t *col = new real_t[4];

	/* Test matrix:
	 *
	 *  [  0  3  6 ]
	 *  [  1  0  7 ]
	 *  [  2  0  8 ]
	 *  [  0  4  9 ]
	 *  [  0  5 10 ]
	 */

	jc[0] = 0; jc[1] = 2; jc[2] = 5; jc[3] = 10;
	ir[0] = 1; ir[1] = 2;
	ir[2] = 0; ir[3] = 3; ir[4] = 4;
	ir[5] = 0; ir[6] = 1; ir[7] = 2; ir[8] = 3; ir[9] = 4;
	for (i = 0; i < 10; i++) val[i] = 1.0+i;

	SparseMatrix A(5, 3, ir, jc, val);

	Indexlist rows(4);

	rows.addNumber(2);
	rows.addNumber(4);
	rows.addNumber(3);
	rows.addNumber(0);

	/* Indexed matrix:
	 *
	 *  [  2  0  8 ]
	 *  [  0  5 10 ]
	 *  [  0  4  9 ]
	 *  [  0  3  6 ]
	 */

	for (j = 0; j < 3; j++)
	{
		fprintf(stdFile, "Column %ld:\n", j);
		A.getCol(j, &rows, 1.0, col);
		for (i = 0; i < 4; i++)
			fprintf(stdFile, " %3.0f\n", col[i]);
	}
				
	delete[] col;
	A.free ();	// or delete[] val,jc,ir;

	return 0;
}

long spGetRow()
{
	long j, i;
	sparse_int_t *ir = new sparse_int_t[10];
	sparse_int_t *jc = new sparse_int_t[6];
	real_t *val = new real_t[10];
	real_t *row = new real_t[4];

	/* Test matrix:
	 *
	 *  [  0  3  4  6  0 ]
	 *  [  1  0  0  7  9 ]
	 *  [  2  0  5  8 10 ]
	 */

	jc[0] = 0; jc[1] = 2; jc[2] = 3; jc[3] = 5; jc[4] = 8; jc[5] = 10;
	ir[0] = 1; ir[1] = 2;
	ir[2] = 0; 
	ir[3] = 0; ir[4] = 2;
	ir[5] = 0; ir[6] = 1; ir[7] = 2; 
	ir[8] = 1; ir[9] = 2;
	for (i = 0; i < 10; i++) val[i] = 1.0+i;

	SparseMatrix A(3, 4, ir, jc, val);

	Indexlist cols(4);

	cols.addNumber(2);
	cols.addNumber(4);
	cols.addNumber(3);
	cols.addNumber(1);

	/* Indexed matrix:
	 *
	 * [  4  0  6  3 ]
	 * [  0  9  7  0 ]
	 * [  5 10  8  0 ]
	 */

	for (j = 0; j < 3; j++)
	{
		A.getRow(j, &cols, 1.0, row);
		for (i = 0; i < 4; i++)
			fprintf(stdFile, " %3.0f", row[i]);
		fprintf(stdFile, "\n");
	}
				
	delete[] row;
	A.free ();

	return 0;
}

long spTimes()
{
	long i;
	sparse_int_t *ir = new sparse_int_t[10];
	sparse_int_t *jc = new sparse_int_t[6];
	real_t *val = new real_t[10];

	real_t *x = new real_t[5*2];
	real_t *y = new real_t[3*2];

	real_t Ax[] = {-23, -11, -26, 42, 74, 99};
	real_t ATy[] = {-63, -69, -222, -423, -359, 272, 126, 663, 1562, 1656};
	real_t err = 0.0;

	for (i = 0; i < 10; i++) x[i] = -4.0+i;

	/* Test matrix:
	 *
	 *  [  0  3  4  6  0 ]
	 *  [  1  0  0  7  9 ]
	 *  [  2  0  5  8 10 ]
	 */

	jc[0] = 0; jc[1] = 2; jc[2] = 3; jc[3] = 5; jc[4] = 8; jc[5] = 10;
	ir[0] = 1; ir[1] = 2;
	ir[2] = 0; 
	ir[3] = 0; ir[4] = 2;
	ir[5] = 0; ir[6] = 1; ir[7] = 2; 
	ir[8] = 1; ir[9] = 2;
	for (i = 0; i < 10; i++) val[i] = 1.0+i;

	SparseMatrix A(3, 5, ir, jc, val);	// reference to ir, jc, val

	A.times(2, 1.0, x, 5, 0.0, y, 3);
	
	for (i = 0; i < 6; i++)
		if (fabs(y[i] - Ax[i]) > err) err = fabs(y[i] - Ax[i]);
	fprintf(stdFile, "Error in sparse A*x: %9.2e\n", err);

	A.transTimes(2, 1.0, y, 3, 0.0, x, 5);

	err = 0.0;
	for (i = 0; i < 10; i++)
		if (fabs(x[i] - ATy[i]) > err) err = fabs(x[i] - ATy[i]);
	fprintf(stdFile, "Error in sparse A'*x: %9.2e\n", err);

	A.free ();	// or delete[] val,ir,jc
	delete[] y;
	delete[] x;

	return 0;
}

long spIndTimes()
{
	const long N = 4;
	long i, j;
	long nRows = 2 * N + 1;
	long nCols = N;
	long nnz = 3 * N;
	sparse_int_t *ir = new sparse_int_t[nnz];
	sparse_int_t *jc = new sparse_int_t[nCols+1];
	real_t *val = new real_t[nnz];
	real_t *xc = new real_t[3*2];
	real_t *yc = new real_t[4*2];
	real_t Ax[] = {0.31, 0.05, 0.06, 0.30, 0.76, 0.20, 0.24, 0.60};
	real_t ATy[] = {0.278, 0.000, 0.548, 0.776, 0.000, 1.208};

	Indexlist rows(4), cols(3), allcols(nCols);

	rows.addNumber(2);
	rows.addNumber(4);
	rows.addNumber(3);
	rows.addNumber(0);

	cols.addNumber(1);
	cols.addNumber(3);
	cols.addNumber(0);

	for (i = 0; i < nCols; i++) 
		allcols.addNumber(i);

	// build test matrix
	for (i = 0; i <= N; i++) jc[i] = 3*i;
	for (j = 0; j < N; j++) 
		for (i = 0; i < 3; i++)
		{
			ir[j*3+i] = 2*j + i;
			val[j*3+i] = 1.0 - 0.1 * (j*3+i);
		}
	SparseMatrix A(nRows, nCols, ir, jc, val);

	fprintf(stdFile, "Test matrix A =\n");
	for (j = 0; j < nRows; j++)
	{
		A.getRow(j, &allcols, 1.0, xc);
		for (i = 0; i < nCols; i++)
			fprintf(stdFile, "%6.2f", xc[i]);
		fprintf(stdFile, "\n");
	}

	for (i = 0; i < 6; i++)
		xc[i] = (1.0 + i) * 0.1;

	A.times(&rows, &cols, 2, 1.0, xc, 3, 0.0, yc, 4, BT_TRUE);

	real_t err = 0.0;
	for (i = 0; i < 8; i++)
		if (fabs(yc[i] - Ax[i]) > err)
			err = fabs(yc[i] - Ax[i]);
	fprintf(stdFile, "Error in sparse indexed A*x: %9.2e\n", err);

	A.transTimes(&rows, &cols, 2, 1.0, yc, 4, 0.0, xc, 3);
	err = 0.0;
	for (i = 0; i < 6; i++)
		if (fabs(xc[i] - ATy[i]) > err)
			err = fabs(xc[i] - ATy[i]);
	fprintf(stdFile, "Error in sparse indexed A'*y: %9.2e\n", err);

	delete[] xc;
	delete[] yc;
	A.free ();

	return 0;
}

long sprGetCol()
{
	long j, i;
	sparse_int_t *ir = new sparse_int_t[10];
	sparse_int_t *jc = new sparse_int_t[4];
	real_t *val = new real_t[10];
	real_t *col = new real_t[4];

	/* Test matrix:
	 *
	 *  [  0  3  6 ]
	 *  [  1  0  7 ]
	 *  [  2  0  8 ]
	 *  [  0  4  9 ]
	 *  [  0  5 10 ]
	 */

	jc[0] = 0; jc[1] = 2; jc[2] = 5; jc[3] = 10;
	ir[0] = 1; ir[1] = 2;
	ir[2] = 0; ir[3] = 3; ir[4] = 4;
	ir[5] = 0; ir[6] = 1; ir[7] = 2; ir[8] = 3; ir[9] = 4;
	for (i = 0; i < 10; i++) val[i] = 1.0+i;

	SparseMatrix Ac(5, 3, ir, jc, val);
	real_t *Acv = Ac.full(); // row major format
	SparseMatrixRow A(5, 3, 3, Acv);
	delete[] Acv;
	Ac.free ();	// or delete[] val,jc,ir;

	Indexlist rows(4);

	rows.addNumber(2);
	rows.addNumber(4);
	rows.addNumber(3);
	rows.addNumber(0);

	/* Indexed matrix:
	 *
	 *  [  2  0  8 ]
	 *  [  0  5 10 ]
	 *  [  0  4  9 ]
	 *  [  0  3  6 ]
	 */

	for (j = 0; j < 3; j++)
	{
		fprintf(stdFile, "Column %ld:\n", j);
		A.getCol(j, &rows, 1.0, col);
		for (i = 0; i < 4; i++)
			fprintf(stdFile, " %3.0f\n", col[i]);
	}
				
	delete[] col;
	A.free ();	// or delete[] val,jc,ir;

	return 0;
}

long sprGetRow()
{
	long j, i;
	sparse_int_t *ir = new sparse_int_t[10];
	sparse_int_t *jc = new sparse_int_t[6];
	real_t *val = new real_t[10];
	real_t *row = new real_t[4];

	/* Test matrix:
	 *
	 *  [  0  3  4  6  0 ]
	 *  [  1  0  0  7  9 ]
	 *  [  2  0  5  8 10 ]
	 */

	jc[0] = 0; jc[1] = 2; jc[2] = 3; jc[3] = 5; jc[4] = 8; jc[5] = 10;
	ir[0] = 1; ir[1] = 2;
	ir[2] = 0; 
	ir[3] = 0; ir[4] = 2;
	ir[5] = 0; ir[6] = 1; ir[7] = 2; 
	ir[8] = 1; ir[9] = 2;
	for (i = 0; i < 10; i++) val[i] = 1.0+i;

	SparseMatrix Ac(3, 5, ir, jc, val);
	real_t *Acv = Ac.full(); // row major format
	SparseMatrixRow A(3, 5, 5, Acv);
	delete[] Acv;
	Ac.free ();	// or delete[] val,jc,ir;

	Indexlist cols(4);

	cols.addNumber(2);
	cols.addNumber(4);
	cols.addNumber(3);
	cols.addNumber(1);

	/* Indexed matrix:
	 *
	 * [  4  0  6  3 ]
	 * [  0  9  7  0 ]
	 * [  5 10  8  0 ]
	 */

	for (j = 0; j < 3; j++)
	{
		A.getRow(j, &cols, 1.0, row);
		for (i = 0; i < 4; i++)
			fprintf(stdFile, " %3.0f", row[i]);
		fprintf(stdFile, "\n");
	}
				
	delete[] row;
	A.free ();

	return 0;
}

long sprTimes()
{
	long i;
	sparse_int_t *ir = new sparse_int_t[10];
	sparse_int_t *jc = new sparse_int_t[6];
	real_t *val = new real_t[10];

	real_t *x = new real_t[5*2];
	real_t *y = new real_t[3*2];

	real_t Ax[] = {-23, -11, -26, 42, 74, 99};
	real_t ATy[] = {-63, -69, -222, -423, -359, 272, 126, 663, 1562, 1656};
	real_t err = 0.0;

	for (i = 0; i < 10; i++) x[i] = -4.0+i;

	/* Test matrix:
	 *
	 *  [  0  3  4  6  0 ]
	 *  [  1  0  0  7  9 ]
	 *  [  2  0  5  8 10 ]
	 */

	jc[0] = 0; jc[1] = 2; jc[2] = 3; jc[3] = 5; jc[4] = 8; jc[5] = 10;
	ir[0] = 1; ir[1] = 2;
	ir[2] = 0; 
	ir[3] = 0; ir[4] = 2;
	ir[5] = 0; ir[6] = 1; ir[7] = 2; 
	ir[8] = 1; ir[9] = 2;
	for (i = 0; i < 10; i++) val[i] = 1.0+i;

	SparseMatrix Ac(3, 5, ir, jc, val);	// reference to ir, jc, val
	real_t *Acv = Ac.full(); // row major format
	SparseMatrixRow A(3, 5, 5, Acv);
	delete[] Acv;
	Ac.free ();	// or delete[] val,jc,ir;

	A.times(2, 1.0, x, 5, 0.0, y, 3);
	
	for (i = 0; i < 6; i++)
		if (fabs(y[i] - Ax[i]) > err) err = fabs(y[i] - Ax[i]);
	fprintf(stdFile, "Error in sparse A*x: %9.2e\n", err);

	A.transTimes(2, 1.0, y, 3, 0.0, x, 5);

	err = 0.0;
	for (i = 0; i < 10; i++)
		if (fabs(x[i] - ATy[i]) > err) err = fabs(x[i] - ATy[i]);
	fprintf(stdFile, "Error in sparse A'*x: %9.2e\n", err);

	A.free ();	// or delete[] val,ir,jc
	delete[] y;
	delete[] x;

	return 0;
}

long sprIndTimes()
{
	const long N = 4;
	long i, j;
	long nRows = 2 * N + 1;
	long nCols = N;
	long nnz = 3 * N;
	sparse_int_t *ir = new sparse_int_t[nnz];
	sparse_int_t *jc = new sparse_int_t[nCols+1];
	real_t *val = new real_t[nnz];
	real_t *xc = new real_t[3*2];
	real_t *yc = new real_t[4*2];
	real_t Ax[] = {0.31, 0.05, 0.06, 0.30, 0.76, 0.20, 0.24, 0.60};
	real_t ATy[] = {0.278, 0.000, 0.548, 0.776, 0.000, 1.208};

	Indexlist rows(4), cols(3), allcols(nCols);

	rows.addNumber(2);
	rows.addNumber(4);
	rows.addNumber(3);
	rows.addNumber(0);

	cols.addNumber(1);
	cols.addNumber(3);
	cols.addNumber(0);

	for (i = 0; i < nCols; i++) 
		allcols.addNumber(i);

	// build test matrix
	for (i = 0; i <= N; i++) jc[i] = 3*i;
	for (j = 0; j < N; j++) 
		for (i = 0; i < 3; i++)
		{
			ir[j*3+i] = 2*j + i;
			val[j*3+i] = 1.0 - 0.1 * (j*3+i);
		}
	SparseMatrix Ac(nRows, nCols, ir, jc, val);
	real_t *Acv = Ac.full(); // row major format
	SparseMatrixRow A(nRows, nCols, nCols, Acv);
	delete[] Acv;
	Ac.free ();	// or delete[] val,jc,ir;

	fprintf(stdFile, "Test matrix A =\n");
	for (j = 0; j < nRows; j++)
	{
		A.getRow(j, &allcols, 1.0, xc);
		for (i = 0; i < nCols; i++)
			fprintf(stdFile, "%6.2f", xc[i]);
		fprintf(stdFile, "\n");
	}

	for (i = 0; i < 6; i++)
		xc[i] = (1.0 + i) * 0.1;

	A.times(&rows, &cols, 2, 1.0, xc, 3, 0.0, yc, 4, BT_TRUE);

	real_t err = 0.0;
	for (i = 0; i < 8; i++)
		if (fabs(yc[i] - Ax[i]) > err)
			err = fabs(yc[i] - Ax[i]);
	fprintf(stdFile, "Error in sparse indexed A*x: %9.2e\n", err);

	A.transTimes(&rows, &cols, 2, 1.0, yc, 4, 0.0, xc, 3);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 2; j++)
			fprintf(stdFile, "%6.2f", ATy[i + j*3]);
		fprintf(stdFile, "\n");
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 2; j++)
			fprintf(stdFile, "%6.2f", xc[i + j*3]);
		fprintf(stdFile, "\n");
	}
	err = 0.0;
	for (i = 0; i < 6; i++)
		if (fabs(xc[i] - ATy[i]) > err)
			err = fabs(xc[i] - ATy[i]);
	fprintf(stdFile, "Error in sparse indexed A'*y: %9.2e\n", err);

	delete[] xc;
	delete[] yc;
	A.free ();

	return 0;
}

int Symmetry()
{
	int i,j;
	real_t *Hv = new real_t[6*6];
	real_t *Z = new real_t[6*3];
	real_t *ZHZd = new real_t[3*3];
	real_t *ZHZs = new real_t[3*3];
	real_t ZHZv[] = {0.144, 0.426, 0.708, 0.426, 1.500, 2.574, 0.708, 2.574, 4.440};
	real_t err;
	SymDenseMat *Hd;
	SymSparseMat *Hs;
	Indexlist *cols = new Indexlist(6);

	for (i = 0; i < 36; i++) Hv[i] = 0.0;
	for (i = 0; i < 6; i++) Hv[i*7] = 1.0 - 0.1 * i;
	for (i = 0; i < 5; i++) Hv[i*7+1] = Hv[i*7+6] = -0.1 * (i+1);

	Hd = new SymDenseMat(6, 6, 6, Hv);	// deep-copy from Hv
	Hs = new SymSparseMat(6, 6, 6, Hv);	// deep-copy from Hv
	Hs->createDiagInfo();

	for (i = 0; i < 6; ++i)
	{
		for (j = 0; j < 6; ++j)
			fprintf (stderr, "%3.3f ", Hv[i*6+j]);
		fprintf (stderr, "\n");
	}
	fprintf (stderr, "\n");
	
	cols->addNumber(3);
	cols->addNumber(0);
	cols->addNumber(4);
	cols->addNumber(1);
	for (i = 0; i < 18; i++) Z[i] = 0.1 * (i+1);

	fprintf (stderr, "\n");
	for (i = 0; i < 6; ++i)
	{
		for (j = 0; j < 3; ++j)
			fprintf (stderr, "%3.3f ", Z[i+j*6]);
		fprintf (stderr, "\n");
	}
	fprintf (stderr, "\n");
	
	Hd->bilinear(cols, 3, Z, 6, ZHZd, 3);

	err = 0.0;
	for (i = 0; i < 9; i++)
		if (fabs(ZHZd[i] - ZHZv[i]) > err)
			err = fabs(ZHZd[i] - ZHZv[i]);
	fprintf(stdFile, "Error in indexed dense bilinear form: %9.2e\n", err);

	Hs->bilinear(cols, 3, Z, 6, ZHZs, 3);
	
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			fprintf (stderr, "%3.3f ", ZHZd[i*3+j]);
		fprintf (stderr, "\n");
	}
	fprintf (stderr, "\n");

	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			fprintf (stderr, "%3.3f ", ZHZv[i*3+j]);
		fprintf (stderr, "\n");
	}
	fprintf (stderr, "\n");

	err = 0.0;
	for (i = 0; i < 9; i++)
		if (fabs(ZHZs[i] - ZHZv[i]) > err)
			err = fabs(ZHZs[i] - ZHZv[i]);
	fprintf(stdFile, "Error in indexed sparse bilinear form: %9.2e\n", err);

	delete cols;	
	delete Hs;
	delete Hd;
	delete[] ZHZs;
	delete[] ZHZd;
	delete[] Z;
	delete[] Hv;

	return 0;
}

/** Test matrix classes */
int main()
{
	hilbert();
	submatrix();
	sumOfSquares();
	indexDenseSubmatrix();

	spGetCol();
	spGetRow();
	spTimes();
	spIndTimes();

	sprGetCol();
	sprGetRow();
	sprTimes();
	sprIndTimes();

	Symmetry();

	return 0;
}


/*
 *	end of file
 */
