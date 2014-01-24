#include <math.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "ChebyshevPolynomials.h"

/* LAPACK generic linear solver, using double precision. */
extern void dgesv_ (int *n, int *nrhs, double *a, int *lda, int *ipiv,
                            double *b, int *ldb, int *info);

void Galerkin_1D_1st_order_ODE(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, j, n;
  int status;
  int gsh[1];
  double **matrix, *matrix_transpose;
  double *rhs;
  double n_d; // so I don't have to keep typing "(double)n"

  int nrhs;
  int *pivot;
  int err;

  nrhs = 1;

  status = CCTK_GroupgshVN(cctkGH, 1, gsh, "spectral_method::x");

  /* So that the matrix of spectral coefficients is square, the number of basis
   * functions by construction must equal the number of grid points used. */

  matrix = (double**) malloc(gsh[0]*sizeof(double*));
  for (i=0; i<gsh[0]; i++)
  {
      matrix[i] = (double*) malloc(gsh[0]*sizeof(double));
  }

  rhs = (double*) malloc(gsh[0]*sizeof(double));
  matrix_transpose = (double*) malloc(gsh[0]*gsh[0]*sizeof(double));
  pivot = (int*) malloc(gsh[0]*sizeof(double));

  /* Do all rows except the boundary condition. */
  for (i=0; i<gsh[0]-1; i++)
  {
    for (n=0; n<gsh[0]; n++)
    {
        n_d = (double)n;

        if (n == 0) {
            matrix[i][n] = 0.0;
        } else {
            matrix[i][n] = n_d * cheby_poly_2(n-1, x[i]);
        }

    }
  }
  for (i=0; i<gsh[0]-1; i++)
  {
      rhs[i] = 2.0;
  }

  /* Boundary condition. */
  for (n = 0; n < gsh[0]; n++) {
    matrix[gsh[0]-1][n] = cheby_poly_1(n, 0.0);
  }
  rhs[gsh[0]-1] = 0.0;


  /* Now transpose the matrix so we can use LAPACK (which is column-major). */
  for (j=0; j<gsh[0]; j++)
  {
      for (i=0; i<gsh[0]; i++)
      {
          matrix_transpose[i + gsh[0]*j] = matrix[i][j];
      }
  }
  dgesv_(&gsh[0], &nrhs, matrix_transpose, &gsh[0], pivot, rhs, &gsh[0], &err);
  if (err != 0)
  {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
              "LAPACK returned non-zero error: %d", err);
  }

  for (i=0; i<gsh[0]; i++)
  {
      sc[i] = rhs[i];
  }

  for (i=0; i<gsh[0]; i++)
  {
      free(matrix[i]);
  }
  free(matrix);
  free(matrix_transpose);
  free(pivot);
  free(rhs);
}
