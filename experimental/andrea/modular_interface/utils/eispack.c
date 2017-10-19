#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "eispack.h"

/******************************************************************************/

int bakvec ( int n, double t[], double e[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    BAKVEC determines eigenvectors by reversing the FIGI transformation.

  Discussion:

    This subroutine forms the eigenvectors of a nonsymmetric tridiagonal
    matrix by back transforming those of the corresponding symmetric
    matrix determined by FIGI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double T[N*3], contains the nonsymmetric matrix.  Its
    subdiagonal is stored in the positions 2:N of the first column,
    its diagonal in positions 1:N of the second column,
    and its superdiagonal in positions 1:N-1 of the third column.
    T(1,1) and T(N,3) are arbitrary.

    Input/output, double E[N].  On input, E(2:N) contains the
    subdiagonal elements of the symmetric matrix.  E(1) is arbitrary.
    On output, the contents of E have been destroyed.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double Z[N*M], contains the eigenvectors.
    On output, they have been transformed as requested.

    Output, int BAKVEC, an error flag.
    0, for normal return,
    2*N+I, if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
    In this case, the symmetric matrix is not similar
    to the original matrix, and the eigenvectors
    cannot be found by this program.
*/
{
  int i;
  int ierr;
  int j;

  ierr = 0;

  if ( m == 0 )
  {
    return ierr;
  }

  e[0] = 1.0;
  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    if ( e[i] == 0.0 )
    {
      if ( t[i+0*3] != 0.0 || t[i-1+2*3] != 0.0 )
      {
        ierr = 2 * n + ( i + 1 );
        return ierr;
      }
      e[i] = 1.0;
    }
    else
    {
      e[i] = e[i-1] * e[i] / t[i-1+2*3];
    }
  }

  for ( j = 0; j < m; j++ )
  {
    for ( i = 1; i < n; i++ )
    {
      z[i+j*n] = z[i+j*n] * e[i];
    }
  }

  return ierr;
}
/******************************************************************************/

void balbak ( int n, int low, int igh, double scale[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    BALBAK determines eigenvectors by undoing the BALANC transformation.

  Discussion:

    This subroutine forms the eigenvectors of a real general matrix by
    back transforming those of the corresponding balanced matrix
    determined by BALANC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2013

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Parlett and Reinsch,
    Numerische Mathematik,
    Volume 13, pages 293-304, 1969.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, column indices determined by BALANC.

    Input, double SCALE[N], contains information determining
    the permutations and scaling factors used by BALANC.

    Input, int M, the number of columns of Z to be
    back-transformed.

    Input/output, double Z[N*M], contains the real and imaginary
    parts of the eigenvectors, which, on return, have been back-transformed.
*/
{
  int i;
  int ii;
  int j;
  int k;
  double s;
  double t;

  if ( m <= 0 )
  {
    return;
  }

  if ( igh != low )
  {
    for ( i = low - 1; i <= igh - 1; i++ )
    {
      for ( j = 0; j < m; j++ )
      {
        z[i+j*n] = z[i+j*n] * scale[i];
      }
    }
  }

  for ( ii = 1; ii <= n; ii++ )
  {
    i = ii;

    if ( i < low || igh < i )
    {
      if ( i < low )
      {
        i = low - ii;
      }

      k = ( int ) ( scale[i-1] );

      if ( k != i )
      {
        for ( j = 0; j < m; j++ )
        {
          t          = z[i-1+j*n];
          z[i-1+j*n] = z[k-1+j*n];
          z[k-1+j*n] = t;
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void bandr ( int n, int mb, double a[], double d[], double e[], double e2[],
  int matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    BANDR reduces a symmetric band matrix to symmetric tridiagonal form.

  Discussion:

    This subroutine reduces a real symmetric band matrix
    to a symmetric tridiagonal matrix using and optionally
    accumulating orthogonal similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int MB, is the (half) band width of the matrix,
    defined as the number of adjacent diagonals, including the principal
    diagonal, required to specify the non-zero portion of the
    lower triangle of the matrix.

    Input/output, double A[N*MB].  On input, contains the lower
    triangle of the symmetric band input matrix stored as an N by MB array.
    Its lowest subdiagonal is stored in the last N+1-MB positions of the first
    column, its next subdiagonal in the last N+2-MB positions of the second
    column, further subdiagonals similarly, and finally its principal diagonal
    in the N positions of the last column.  Contents of storages not part of
    the matrix are arbitrary.  On output, A has been destroyed, except for
    its last two columns which contain a copy of the tridiagonal matrix.

    Output, double D[N], the diagonal elements of the tridiagonal
    matrix.

    Output, double E[N], the subdiagonal elements of the tridiagonal
    matrix in E(2:N).  E(1) is set to zero.

    Output, double E2[N], contains the squares of the corresponding
    elements of E.  E2 may coincide with E if the squares are not needed.

    Input, logical MATZ, should be set to TRUE if the transformation matrix is
    to be accumulated, and to FALSE otherwise.

    Output, double Z[N*N], the orthogonal transformation matrix
    produced in the reduction if MATZ has been set to TRUE.  Otherwise, Z is
    not referenced.
*/
{
  double b1;
  double b2;
  double c2;
  double dmin;
  double dminrt;
  double f1;
  double f2;
  double g;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  int jj;
  int k;
  int kr;
  int l;
  int m1;
  int maxl;
  int maxr;
  int mr;
  int r;
  int r1;
  double s2;
  double u;
  int ugl;

  dmin = r8_epsilon ( );
  dminrt = sqrt ( dmin );
/*
  Initialize the diagonal scaling matrix.
*/
  for ( i = 0; i < n; i++ )
  {
    d[i] = 1.0;
  }

  if ( matz )
  {
    r8mat_identity ( n, z );
  }
/*
  Is input matrix diagonal?
*/
  if ( mb == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      d[i] = a[i+(mb-1)*n];
      e[i] = 0.0;
      e2[i] = 0.0;
    }
    return;
  }

  m1 = mb - 1;

  if ( m1 != 1 )
  {
    for ( k = 1; k <= n - 2; k++ )
    {
      maxr = i4_min ( m1, n - k );
      for ( r1 = 2; r1 <= maxr; r1++ )
      {
        r = maxr + 2 - r1;
        kr = k + r;
        mr = mb - r;
        g = a[kr-1+(mr-1)*n];
        a[kr-2+0*n] = a[kr-2+mr*n];
        ugl = k;

        for ( j = kr; j <= n; j = j + m1 )
        {
          j1 = j - 1;
          j2 = j1 - 1;

          if ( g == 0.0 )
          {
            break;
          }

          b1 = a[j1-1+0*n] / g;
          b2 = b1 * d[j1-1] / d[j-1];
          s2 = 1.0 / ( 1.0 + b1 * b2 );

          if ( s2 < 0.5 )
          {
            b1 = g / a[j1-1+0*n];
            b2 = b1 * d[j-1] / d[j1-1];
            c2 = 1.0 - s2;
            d[j1-1] = c2 * d[j1-1];
            d[j-1] = c2 * d[j-1];
            f1 = 2.0 * a[j-1+(m1-1)*n];
            f2 = b1 * a[j1-1+(mb-1)*n];
            a[j-1+(m1-1)*n] = -b2 * ( b1 * a[j-1+(m1-1)*n] - a[j-1+(mb-1)*n] ) - f2 + a[j-1+(m1-1)*n];
            a[j1-1+(mb-1)*n] = b2 * ( b2 * a[j-1+(mb-1)*n] + f1 ) + a[j1-1+(mb-1)*n];
            a[j-1+(mb-1)*n] = b1 * ( f2 - f1 ) + a[j-1+(mb-1)*n];

            for ( l = ugl; l <= j2; l++ )
            {
              i2 = mb - j + l;
              u = a[j1-1+i2*n] + b2 * a[j-1+(i2-1)*n];
              a[j-1+(i2-1)*n] = -b1 * a[j1-1+i2*n] + a[j-1+(i2-1)*n];
              a[j1-1+i2*n] = u;
            }

            ugl = j;
            a[j1-1+0*n] = a[j1+0*n] + b2 * g;

            if ( j != n )
            {
              maxl = i4_min ( m1, n - j1 );

              for ( l = 2; l <= maxl; l++ )
              {
                i1 = j1 + l;
                i2 = mb - l;
                u = a[i1-1+(i2-1)*n] + b2 * a[i1-1+i2*n];
                a[i1-1+i2*n] = -b1 * a[i1-1+(i2-1)*n] + a[i1-1+i2*n];
                a[i1-1+(i2-1)*n] = u;
              }

              i1 = j + m1;

              if ( i1 <= n )
              {
                g = b2 * a[i1-1+0*n];
              }
            }

            if ( matz )
            {
              for ( l = 1; l <= n; l++ )
              {
                u = z[l-1+(j1-1)*n] + b2 * z[l-1+(j-1)*n];
                z[l-1+(j-1)*n] = -b1 * z[l-1+(j1-1)*n] + z[l-1+(j-1)*n];
                z[l-1+(j1-1)*n] = u;
              }
            }
          }
          else
          {
            u = d[j1-1];
            d[j1-1] = s2 * d[j-1];
            d[j-1] = s2 * u;
            f1 = 2.0 * a[j-1+(m1-1)*n];
            f2 = b1 * a[j-1+(mb-1)*n];
            u = b1 * ( f2 - f1 ) + a[j1-1+(mb-1)*n];
            a[j-1+(m1-1)*n] = b2 * ( b1 * a[j-1+(m1-1)*n] - a[j1-1+(mb-1)*n] ) + f2 - a[j-1+(m1-1)*n];
            a[j1-1+(mb-1)*n] = b2 * ( b2 * a[j1-1+(mb-1)*n] + f1 ) + a[j-1+(mb-1)*n];
            a[j-1+(mb-1)*n] = u;

            for ( l = ugl; l <= j2; l++ )
            {
              i2 = mb - j + l;
              u = b2 * a[j1-1+i2*n] + a[j-1+(i2-1)*n];
              a[j-1+(i2-1)*n] = -a[j1-1+i2*n] + b1 * a[j-1+(i2-1)*n];
              a[j1-1+i2*n] = u;
            }

            ugl = j;
            a[j1-1+0*n] = b2 * a[j1-1+0*n] + g;

            if ( j != n )
            {
              maxl = i4_min ( m1, n - j1 );

              for ( l = 2; l <= maxl; l++ )
              {
                i1 = j1 + l;
                i2 = mb - l;
                u = b2 * a[i1-1+(i2-1)*n] + a[i1-1+i2*n];
                a[i1-1+i2*n] = -a[i1-1+(i2-1)*n] + b1 * a[i1-1+i2*n];
                a[i1-1+(i2-1)*n] = u;
              }

              i1 = j + m1;

              if ( i1 <= n )
              {
                g = a[i1-1+0*n];
                a[i1-1+0*n] = b1 * a[i1-1+0*n];
              }
            }

            if ( matz )
            {
              for ( l = 1; l <= n; l++ )
              {
                u = b2 * z[l-1+(j1-1)*n] + z[l-1+(j-1)*n];
                z[l-1+(j-1)*n] = -z[l-1+(j1-1)*n] + b1 * z[l-1+(j-1)*n];
                z[l-1+(j1-1)*n] = u;
              }
            }
          }
        }
      }
/*
  Rescale to avoid underflow or overflow.
*/
      if ( ( k % 64 ) == 0 )
      {
        for ( j = k; j <= n; j++ )
        {
          if ( d[j-1] < dmin )
          {
            maxl = i4_max ( 1, mb + 1 - j );

            for ( jj = maxl; jj <= m1; jj++ )
            {
              a[j-1+(jj-1)*n] = dminrt * a[j-1+(jj-1)*n];
            }

            if ( j != n )
            {
              maxl = i4_min ( m1, n - j );

              for ( l = 1; l <= maxl; l++ )
              {
                i1 = j + l;
                i2 = mb - l;
                a[i1-1+(i2-1)*n] = dminrt * a[i1-1+(i2-1)*n];
              }
            }

            if ( matz )
            {
              for ( i = 1; i <= n; i++ )
              {
                z[i-1+(j-1)*n] = dminrt * z[i-1+(j-1)*n];
              }
            }

            a[j-1+(mb-1)*n] = dmin * a[j-1+(mb-1)*n];
            d[j-1] = d[j-1] / dmin;
          }
        }
      }
    }
  }
/*
  Form square root of scaling matrix.
*/
  for ( i = 1; i < n; i++ )
  {
    e[i] = sqrt ( d[i] );
  }
  if ( matz )
  {
    for ( j = 1; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i+j*n] = z[i+j*n] * e[j];
      }
    }
  }

  u = 1.0;

  for ( j = 1; j < n; j++ )
  {
    a[j+(m1-1)*n] = u * e[j] * a[j+(m1-1)*n];
    u = e[j];
    e2[j] = a[j+(m1-1)*n] * a[j+(m1-1)*n];
    a[j+(mb-1)*n] = d[j] * a[j+(mb-1)*n];
    d[j] = a[j+(mb-1)*n];
    e[j] = a[j+(m1-1)*n];
  }

  d[0] = a[0+(mb-1)*n];
  e[0] = 0.0;
  e2[0] = 0.0;

  return;
}
/******************************************************************************/

void cbabk2 ( int n, int low, int igh, double scale[], int m, double zr[],
  double zi[] )

/******************************************************************************/
/*
  Purpose:

    CBABK2 finds eigenvectors by undoing the CBAL transformation.

  Discussion:

    This subroutine forms the eigenvectors of a complex general
    matrix by back transforming those of the corresponding
    balanced matrix determined by CBAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, values determined by CBAL.

    Input, double SCALE[N], information determining the permutations
    and scaling factors used by CBAL.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double ZR[N*M], ZI[N*M].  On input, the real
    and imaginary parts, respectively, of the eigenvectors to be back
    transformed in their first M columns.  On output, the transformed
    eigenvectors.
*/
{
  int i;
  int ii;
  int j;
  int k;
  double s;

  if ( m == 0 )
  {
    return;
  }

  if ( igh != low )
  {
    for ( i = low; i <= igh; i++ )
    {
      s = scale[i];
      for ( j = 0; j < m; j++ )
      {
        zr[i+j*n] = zr[i+j*n] * s;
        zi[i+j*n] = zi[i+j*n] * s;
      }
    }
  }

  for ( ii = 0; ii < n; ii++ )
  {
    i = ii;
    if ( i < low || igh < i )
    {
      if ( i < low )
      {
        i = low - ii;
      }

      k = scale[i];

      if ( k != i )
      {
        for ( j = 0; j < m; j++ )
        {
          s         = zr[i+j*n];
          zr[i+j*n] = zr[k+j*n];
          zr[k+j*n] = s;
          s         = zi[i+j*n];
          zi[i+j*n] = zi[k+j*n];
          zi[k+j*n] = s;
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void csroot ( double xr, double xi, double *yr, double *yi )

/******************************************************************************/
/*
  Purpose:

    CSROOT computes the complex square root of a complex quantity.

  Discussion:

    The branch of the square function is chosen so that
      0.0 <= YR
    and
      sign ( YI ) == sign ( XI )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, double XR, XI, the real and imaginary parts of the
    quantity whose square root is desired.

    Output, double *YR, *YI, the real and imaginary parts of the
    square root.
*/
{
  double s;
  double ti;
  double tr;

  tr = xr;
  ti = xi;
  s = sqrt ( 0.5 * ( pythag ( tr, ti ) + r8_abs ( tr ) ) );

  if ( 0.0 <= tr )
  {
    *yr = s;
  }

  if ( ti < 0.0 )
  {
    s = -s;
  }

  if ( tr <= 0.0 )
  {
    *yi = s;
  }

  if ( tr < 0.0 )
  {
    *yr = 0.5 * ( ti / *yi );
  }
  else if ( 0.0 < tr )
  {
    *yi = 0.5 * ( ti / *yr );
  }
  return;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

double pythag ( double a, double b )

/******************************************************************************/
/*
  Purpose:

    PYTHAG computes SQRT ( A * A + B * B ) carefully.

  Discussion:

    The formula

      PYTHAG = sqrt ( A * A + B * B )

    is reasonably accurate, but can fail if, for example, A^2 is larger
    than the machine overflow.  The formula can lose most of its accuracy
    if the sum of the squares is very large or very small.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Modified:

    08 November 2012

  Parameters:

    Input, double A, B, the two legs of a right triangle.

    Output, double PYTHAG, the length of the hypotenuse.
*/
{
  double p;
  double r;
  double s;
  double t;
  double u;

  p = r8_max ( r8_abs ( a ), r8_abs ( b ) );

  if ( p != 0.0 )
  {
    r = r8_min ( r8_abs ( a ), r8_abs ( b ) ) / p;
    r = r * r;

    while ( 1 )
    {
      t = 4.0 + r;

      if ( t == 4.0 )
      {
        break;
      }

      s = r / t;
      u = 1.0 + 2.0 * s;
      p = u * p;
      r = ( s / u ) * ( s / u ) * r;
    }
  }
  return p;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
}
/******************************************************************************/

void r8mat_identity  ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY sets an R8MAT to the identity matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double A[N*N], the N by N identity matrix.
*/
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM_NEW multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

int rs ( int n, double a[], double w[], int matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RS computes eigenvalues and eigenvectors of real symmetric matrix.

  Discussion:

    This subroutine calls the recommended sequence of
    subroutines from the eigensystem subroutine package (eispack)
    to find the eigenvalues and eigenvectors (if desired)
    of a real symmetric matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the real symmetric matrix.

    Input, int MATZ, is zero if only eigenvalues are desired,
    and nonzero if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N*N], contains the eigenvectors, if MATZ
    is nonzero.

    Output, int RS, is set equal to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  if ( matz == 0 )
  {
    fv1 = ( double * ) malloc ( n * sizeof ( double ) );
    fv2 = ( double * ) malloc ( n * sizeof ( double ) );

    tred1 ( n, a, w, fv1, fv2 );
    ierr = tqlrat ( n, w, fv2 );

    free ( fv1 );
    free ( fv2 );
  }
  else
  {
    fv1 = ( double * ) malloc ( n * sizeof ( double ) );

    tred2 ( n, a, w, fv1, z );
    ierr = tql2 ( n, w, fv1, z );

    free ( fv1 );
  }

  return ierr;
}
/******************************************************************************/

int rsb ( int n, int mb, double a[], double w[], int matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.

  Discussion:

    This subroutine calls the recommended sequence of
    subroutines from the eigensystem subroutine package (eispack)
    to find the eigenvalues and eigenvectors (if desired)
    of a real symmetric band matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int MB, the half band width of the matrix,
    defined as the number of adjacent diagonals, including the principal
    diagonal, required to specify the non-zero portion of the lower triangle
    of the matrix.

    Input, double A[N*MB], contains the lower triangle of the real
    symmetric band matrix.  Its lowest subdiagonal is stored in the last N+1-MB
    positions of the first column, its next subdiagonal in the last
    N+2-MB positions of the second column, further subdiagonals similarly,
    and finally its principal diagonal in the N positions of the last
    column.  Contents of storages not part of the matrix are arbitrary.

    Input, int MATZ, is zero if only eigenvalues are desired,
    and nonzero if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N*N], contains the eigenvectors, if MATZ
    is nonzero.

    Output, int BANDR, is set to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;
  int tf;

  if ( mb <= 0 )
  {
    ierr = 12 * n;
    return ierr;
  }

  if ( n < mb )
  {
    ierr = 12 * n;
    return ierr;
  }

  if ( matz == 0 )
  {
    fv1 = ( double * ) malloc ( n * sizeof ( double ) );
    fv2 = ( double * ) malloc ( n * sizeof ( double ) );
    tf = 0;

    bandr ( n, mb, a, w, fv1, fv2, tf, z );

    ierr = tqlrat ( n, w, fv2 );

    free ( fv1 );
    free ( fv2 );
  }
  else
  {
    fv1 = ( double * ) malloc ( n * sizeof ( double ) );
    tf = 1;

    bandr ( n, mb, a, w, fv1, fv1, tf, z );

    ierr = tql2 ( n, w, fv1, z );

    free ( fv1 );
  }
  return ierr;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

int tql2 ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.

  Discussion:

    This subroutine finds the eigenvalues and eigenvectors of a symmetric
    tridiagonal matrix by the QL method.  The eigenvectors of a full
    symmetric matrix can also be found if TRED2 has been used to reduce this
    full matrix to tridiagonal form.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Bowdler, Martin, Reinsch, Wilkinson,
    TQL2,
    Numerische Mathematik,
    Volume 11, pages 293-306, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D[N].  On input, the diagonal elements of
    the matrix.  On output, the eigenvalues in ascending order.  If an error
    exit is made, the eigenvalues are correct but unordered for indices
    1,2,...,IERR-1.

    Input/output, double E[N].  On input, E(2:N) contains the
    subdiagonal elements of the input matrix, and E(1) is arbitrary.
    On output, E has been destroyed.

    Input, double Z[N*N].  On input, the transformation matrix
    produced in the reduction by TRED2, if performed.  If the eigenvectors of
    the tridiagonal matrix are desired, Z must contain the identity matrix.
    On output, Z contains the orthonormal eigenvectors of the symmetric
    tridiagonal (or full) matrix.  If an error exit is made, Z contains
    the eigenvectors associated with the stored eigenvalues.

    Output, int TQL2, error flag.
    0, normal return,
    J, if the J-th eigenvalue has not been determined after
    30 iterations.
*/
{
  double c;
  double c2;
  double c3;
  double dl1;
  double el1;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int j;
  int k;
  int l;
  int l1;
  int l2;
  int m;
  int mml;
  double p;
  double r;
  double s;
  double s2;
  double t;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }

  f = 0.0;
  tst1 = 0.0;
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    j = 0;
    h = r8_abs ( d[l] ) + r8_abs ( e[l] );
    tst1 = r8_max ( tst1, h );
/*
  Look for a small sub-diagonal element.
*/
    for ( m = l; m < n; m++ )
    {
      tst2 = tst1 + r8_abs ( e[m] );
      if ( tst2 == tst1 )
      {
        break;
      }
    }

    if ( m != l )
    {
      for ( ; ; )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
/*
  Form shift.
*/
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * e[l] );
        r = pythag ( p, 1.0 );
        d[l] = e[l] / ( p + r8_sign ( p ) * r8_abs ( r ) );
        d[l1] = e[l] * ( p + r8_sign ( p ) * r8_abs ( r ) );
        dl1 = d[l1];
        h = g - d[l];
        for ( i = l2; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
/*
  QL transformation.
*/
        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          i = m - ii;
          g = c * e[i];
          h = c * p;
          r = pythag ( p, e[i] );
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );
/*
  Form vector.
*/
          for ( k = 0; k < n; k++ )
          {
            h = z[k+(i+1)*n];
            z[k+(i+1)*n] = s * z[k+i*n] + c * h;
            z[k+i*n] = c * z[k+i*n] - s * h;
          }
        }
        p = - s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + r8_abs ( e[l] );

        if ( tst2 <= tst1 )
        {
          break;
        }
      }
    }
    d[l] = d[l] + f;
  }
/*
  Order eigenvalues and eigenvectors.
*/
  for ( ii = 1; ii < n; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i];
    for ( j = ii; j < n; j++ )
    {
      if ( d[j] < p )
      {
        k = j;
        p = d[j];
      }
    }

    if ( k != i )
    {
      d[k] = d[i];
      d[i] = p;
      for ( j = 0; j < n; j++ )
      {
        t        = z[j+i*n];
        z[j+i*n] = z[j+k*n];
        z[j+k*n] = t;
      }
    }
  }
  return ierr;
}
/******************************************************************************/

int tqlrat ( int n, double d[], double e2[] )

/******************************************************************************/
/*
  Purpose:

    TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    This subroutine finds the eigenvalues of a symmetric
    tridiagonal matrix by the rational QL method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Christian Reinsch,
    Algorithm 464, TQLRAT,
    Communications of the ACM,
    Volume 16, page 689, 1973.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D[N].  On input, D contains the diagonal
    elements of the matrix.  On output, D contains the eigenvalues in ascending
    order.  If an error exit was made, then the eigenvalues are correct
    in positions 1 through IERR-1, but may not be the smallest eigenvalues.

    Input/output, double E2[N], contains in positions 2 through N
    the squares of the subdiagonal elements of the matrix.  E2(1) is
    arbitrary.  On output, E2 has been overwritten by workspace
    information.

    Output, int TQLRAT, error flag.
    0, for no error,
    J, if the J-th eigenvalue could not be determined after 30 iterations.
*/
{
  double b;
  double c;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int j;
  int l;
  int l1;
  int m;
  int mml;
  double p;
  double r;
  double s;
  double t;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e2[i-1] = e2[i];
  }

  f = 0.0;
  t = 0.0;
  e2[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
     j = 0;
     h = r8_abs ( d[l] ) + sqrt ( e2[l] );

     if ( t <= h )
     {
       t = h;
       b = r8_abs ( t ) * r8_epsilon ( );
       c = b * b;
     }
/*
  Look for small squared sub-diagonal element.
*/
    for ( m = l; m < n; m++ )
    {
      if ( e2[m] <= c )
      {
        break;
      }
    }

    if ( m != l )
    {
      for ( ; ; )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
/*
  Form shift.
*/
        l1 = l + 1;
        s = sqrt ( e2[l] );
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * s );
        r = pythag ( p, 1.0 );
        d[l] = s / ( p + r8_abs ( r ) * r8_sign ( p ) );
        h = g - d[l];
        for ( i = l1; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
/*
  Rational QL transformation.
*/
        g = d[m];
        if ( g == 0.0 )
        {
          g = b;
        }

        h = g;
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          i = m - ii;
          p = g * h;
          r = p + e2[i];
          e2[i+1] = s * r;
          s = e2[i] / r;
          d[i+1] = h + s * ( h + d[i] );
          g = d[i] - e2[i] / g;
          if ( g == 0.0 )
          {
            g = b;
          }
          h = g * p / r;
        }
        e2[l] = s * g;
        d[l] = h;
/*
  Guard against underflow in convergence test.
*/
        if ( h == 0.0 )
        {
          break;
        }

        if ( r8_abs ( e2[l] ) <= r8_abs ( c / h ) )
        {
          break;
        }

        e2[l] = h * e2[l];

        if ( e2[l] == 0.0 )
        {
          break;
        }
      }
    }

    p = d[l] + f;
/*
  Order the eigenvalues.
*/
    for ( i = l; 0 <= i; i-- )
    {
      if ( i == 0 )
      {
        d[i] = p;
        break;
      }
      else if ( d[i-1] <= p )
      {
        d[i] = p;
        break;
      }
      d[i] = d[i-1];
    }
  }

  return ierr;
}
/******************************************************************************/

void tred1 ( int n, double a[], double d[], double e[], double e2[] )

/******************************************************************************/
/*
  Purpose:

    TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.

  Discussion:

    The routine reduces a real symmetric matrix to a symmetric
    tridiagonal matrix using orthogonal similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, Reinsch, Wilkinson,
    TRED1,
    Numerische Mathematik,
    Volume 11, pages 181-195, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix A.

    Input/output, double A[N*N], on input, contains the real
    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
    On output, A contains information about the orthogonal transformations
    used in the reduction in its strict lower triangle.
    The full upper triangle of A is unaltered.

    Output, double D[N], contains the diagonal elements of the
    tridiagonal matrix.

    Output, double E[N], contains the subdiagonal elements of the
    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.

    Output, double E2[N], contains the squares of the corresponding
    elements of E.  E2 may coincide with E if the squares are not needed.
*/
{
  double f;
  double g;
  double h;
  int i;
  int ii;
  int j;
  int k;
  int l;
  double scale;

  for ( j = 0; j < n; j++ )
  {
    d[j] = a[n-1+j*n];
  }

  for ( i = 0; i < n; i++ )
  {
    a[n-1+i*n] = a[i+i*n];
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    l = i - 1;
    h = 0.0;
/*
  Scale row.
*/
    scale = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      scale = scale + r8_abs ( d[k] );
    }

    if ( scale == 0.0 )
    {
      for ( j = 0; j <= l; j++ )
      {
        d[j]     = a[l+j*n];
        a[l+j*n] = a[i+j*n];
        a[i+j*n] = 0.0;
      }

      e[i] = 0.0;
      e2[i] = 0.0;
      continue;
    }

    for ( k = 0; k <= l; k++ )
    {
      d[k] = d[k] / scale;
    }

    for ( k = 0; k <= l; k++ )
    {
      h = h + d[k] * d[k];
    }

    e2[i] = h * scale * scale;
    f = d[l];
    g = - sqrt ( h ) * r8_sign ( f );
    e[i] = scale * g;
    h = h - f * g;
    d[l] = f - g;

    if ( 0 <= l )
    {
/*
  Form A * U.
*/
      for ( k = 0; k <= l; k++ )
      {
        e[k] = 0.0;
      }

      for ( j = 0; j <= l; j++ )
      {
        f = d[j];
        g = e[j] + a[j+j*n] * f;

        for ( k = j + 1; k <= l; k++ )
        {
          g = g + a[k+j*n] * d[k];
          e[k] = e[k] + a[k+j*n] * f;
        }
        e[j] = g;
      }
/*
  Form P.
*/
      f = 0.0;
      for ( j = 0; j <= l; j++ )
      {
        e[j] = e[j] / h;
        f = f + e[j] * d[j];
      }

      h = f / ( h + h );
/*
  Form Q.
*/
      for ( j = 0; j <= l; j++ )
      {
        e[j] = e[j] - h * d[j];
      }
/*
  Form reduced A.
*/
      for ( j = 0; j <= l; j++ )
      {
        f = d[j];
        g = e[j];
        for ( k = j; k <= l; k++ )
        {
          a[k+j*n] = a[k+j*n] - f * e[k] - g * d[k];
        }
      }
    }

    for ( j = 0; j <= l; j++ )
    {
      f        = d[j];
      d[j]     = a[l+j*n];
      a[l+j*n] = a[i+j*n];
      a[i+j*n] = f * scale;
    }
  }
  return;
}
/******************************************************************************/

void tred2 ( int n, double a[], double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.

  Discussion:

    This subroutine reduces a real symmetric matrix to a
    symmetric tridiagonal matrix using and accumulating
    orthogonal similarity transformations.

    A and Z may coincide, in which case a single storage area is used
    for the input of A and the output of Z.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, Reinsch, Wilkinson,
    TRED2,
    Numerische Mathematik,
    Volume 11, pages 181-195, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the real symmetric input matrix.  Only the
    lower triangle of the matrix need be supplied.

    Output, double D[N], the diagonal elements of the tridiagonal
    matrix.

    Output, double E[N], contains the subdiagonal elements of the
    tridiagonal matrix in E(2:N).  E(1) is set to zero.

    Output, double Z[N*N], the orthogonal transformation matrix
    produced in the reduction.
*/
{
  double f;
  double g;
  double h;
  double hh;
  int i;
  int ii;
  int j;
  int k;
  int l;
  double scale;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j; i < n; i++ )
    {
      z[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = a[n-1+j*n];
  }

  for ( i = n - 1; 1 <= i; i-- )
  {
    l = i - 1;
    h = 0.0;
/*
  Scale row.
*/
    scale = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      scale = scale + r8_abs ( d[k] );
    }

    if ( scale == 0.0 )
    {
      e[i] = d[l];

      for ( j = 0; j <= l; j++ )
      {
        d[j]     = z[l+j*n];
        z[i+j*n] = 0.0;
        z[j+i*n] = 0.0;
      }
      d[i] = 0.0;
      continue;
    }

    for ( k = 0; k <= l; k++ )
    {
      d[k] = d[k] / scale;
    }

    h = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      h = h + d[k] * d[k];
    }

    f = d[l];
    g = - sqrt ( h ) * r8_sign ( f );
    e[i] = scale * g;
    h = h - f * g;
    d[l] = f - g;
/*
  Form A*U.
*/
    for ( k = 0; k <= l; k++ )
    {
      e[k] = 0.0;
    }

    for ( j = 0; j <= l; j++ )
    {
      f = d[j];
      z[j+i*n] = f;
      g = e[j] + z[j+j*n] * f;

      for ( k = j + 1; k <= l; k++ )
      {
        g = g + z[k+j*n] * d[k];
        e[k] = e[k] + z[k+j*n] * f;
      }
      e[j] = g;
    }
/*
  Form P.
*/
    for ( k = 0; k <= l; k++ )
    {
      e[k] = e[k] / h;
    }
    f = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      f = f + e[k] * d[k];
    }
    hh = 0.5 * f / h;
/*
  Form Q.
*/
    for ( k = 0; k <= l; k++ )
    {
      e[k] = e[k] - hh * d[k];
    }
/*
  Form reduced A.
*/
    for ( j = 0; j <= l; j++ )
    {
      f = d[j];
      g = e[j];

      for ( k = j; k <= l; k++ )
      {
        z[k+j*n] = z[k+j*n] - f * e[k] - g * d[k];
      }
      d[j] = z[l+j*n];
      z[i+j*n] = 0.0;
    }
    d[i] = h;
  }
/*
  Accumulation of transformation matrices.
*/
  for ( i = 1; i < n; i++ )
  {
    l = i - 1;
    z[n-1+l*n] = z[l+l*n];
    z[l+l*n] = 1.0;
    h = d[i];

    if ( h != 0.0 )
    {
      for ( k = 0; k <= l; k++ )
      {
        d[k] = z[k+i*n] / h;
      }
      for ( j = 0; j <= l; j++ )
      {
        g = 0.0;
        for ( k = 0; k <= l; k++ )
        {
          g = g + z[k+i*n] * z[k+j*n];
        }
        for ( k = 0; k <= l; k++ )
        {
          z[k+j*n] = z[k+j*n] - g * d[k];
        }
      }
    }
    for ( k = 0; k <= l; k++ )
    {
      z[k+i*n] = 0.0;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = z[n-1+j*n];
  }

  for ( j = 0; j < n - 1; j++ )
  {
    z[n-1+j*n] = 0.0;
  }
  z[n-1+(n-1)*n] = 1.0;

  e[0] = 0.0;

  return;
}
