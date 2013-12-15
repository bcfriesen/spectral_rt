#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <error_codes.h>
#include <chebyshev_polynomials.h>
#include <approx_result.h>

extern void dgesv_ (int *n, int *nrhs, double *a, int *lda, int *ipiv,
                    double *b, int *ldb, int *info);

const int SIZE = 10;
const int SIZE_HIGH_RES = 100;

int
main ()
{

  int i, j, k;
  int err;
  double rhs[SIZE];
  const double lower_limit = -1.0;
  const double upper_limit = +1.0;
  double integration_points[SIZE];
  double integration_weights[SIZE];
  double matrix[SIZE][SIZE], matrix_transpose[SIZE * SIZE];
  double total;
  int c1, c2;
  int pivot[SIZE];
  double approx_solution[SIZE_HIGH_RES];
  double high_res_ordinate[SIZE_HIGH_RES];
  double exact_solution[SIZE_HIGH_RES];
  gsl_integration_glfixed_table *integration_table =
    gsl_integration_glfixed_table_alloc (SIZE);
  FILE *fp;
  const char solutions_filename[] = "solutions.dat";
  const char basis_functions_filename[] = "basis_functions.dat";
  const char coeffs_filename[] = "coefficients.dat";

  /* Sanity checks for grid sizes. */
  assert (SIZE > 0);
  assert (SIZE_HIGH_RES > 0);

  /* Calculate integration points and weights using Gauss-Legendre rules. */
  for (i = 0; i < SIZE; ++i)
    {
      err =
        gsl_integration_glfixed_point (lower_limit, upper_limit, i,
                                       &integration_points[i],
                                       &integration_weights[i],
                                       integration_table);
      if (err != GSL_SUCCESS)
        {
          printf
            ("ERROR: could not calculate Gauss-Legendre integration points and weights!\n");
          return (-1);
        }
    }
  gsl_integration_glfixed_table_free (integration_table);


  /* Set initial conditions in both the matrix and the RHS vector. */
  j = SIZE - 1;
  rhs[j] = 0.0;
  for (i = 0; i < SIZE; ++i)
    {
      matrix[j][i] = chebyshev_poly_1 (i, 0.0);
    }


  /* Fill in the rest of the matrix elements A_ji. */
  for (j = 0; j < SIZE - 1; ++j)
    {
      /* Start the column index "i" at 1 instead of 0 because each matrix element contains U_{i-1}. */
      for (i = 1; i < SIZE; ++i)
        {
          total = 0.0;
          for (k = 0; k < SIZE; ++k)
            {
              total =
                total + integration_weights[k] * chebyshev_poly_1 (j,
                                                                   integration_points
                                                                   [k]) *
                (double) i *chebyshev_poly_2 (i - 1, integration_points[k]);
            }
          matrix[j][i] = total;
        }
    }

  /* Fill in the rest of the RHS vector. */
  for (j = 0; j < SIZE - 1; ++j)
    {
      total = 0.0;
      for (k = 0; k < SIZE; ++k)
        {
          total =
            total + integration_weights[k] * chebyshev_poly_1 (j,
                                                               integration_points
                                                               [k]) *
            integration_points[k];
        }
      rhs[j] = 2.0 * total;
    }

  /* Now transpose the matrix so we can use LAPACK routines. */
  for (j = 0; j < SIZE; ++j)
    {
      for (i = 0; i < SIZE; ++i)
        {
          matrix_transpose[i + SIZE * j] = matrix[i][j];
        }
    }

  /* Now solve the linear system for the coefficients. The input vector is overwritten with the results. */
  c1 = SIZE;
  c2 = 1;
  dgesv_ (&c1, &c2, matrix_transpose, &c1, pivot, rhs, &c1, &err);
  if (err != SUCCESS)
    {
      printf
        ("ERROR: LAPACK routine dgesv() failed with the error code: %d\n",
         err);
      return -1;
    }

  /* Now that we have the coefficients, let's see how the numerical solution
   * compares at points other than just the points used above. */
  for (i = 0; i < SIZE_HIGH_RES; ++i)
    {
      high_res_ordinate[i] =
        lower_limit + (double) i *(upper_limit -
                                   lower_limit) / (double) (SIZE_HIGH_RES -
                                                            1);
    }

  for (i = 0; i < SIZE_HIGH_RES; ++i)
    {
      exact_solution[i] = pow (high_res_ordinate[i], 2);
      approx_solution[i] = approx_result (high_res_ordinate[i], rhs, SIZE);
    }

  /* Write the exact and numerical solutions to a file. */
  /* TODO: make the file name an adjustable parameter. */

  fp = fopen (solutions_filename, "w");
  if (fp == NULL)
    {
      fprintf (stderr, "Can't open file '%s' for writing!\n",
               solutions_filename);
      return -1;
    }
  fprintf (fp, "%3s %20s %20s %20s\n", "#", "ORDINATE", "EXACT SOLUTION",
           "NUMERICAL SOLUTION");
  fprintf (fp, "%3s %20s %20s %20s\n", "#", "--------", "--------------",
           "------------------");
  for (i = 0; i < SIZE_HIGH_RES; ++i)
    {
      fprintf (fp, "%3s %20.4e %20.4e %20.4e\n", "   ", high_res_ordinate[i],
               exact_solution[i], approx_solution[i]);
    }
  fclose (fp);
  printf ("Saved numerical solution to file '%s'.\n", solutions_filename);

  /* Write the basis functions to a file. */
  /* TODO: make the file name an adjustable parameter. */

  fp = fopen (basis_functions_filename, "w");
  if (fp == NULL)
    {
      fprintf (stderr, "Can't open file '%s' for writing!\n",
               basis_functions_filename);
      return -1;
    }

  fprintf (fp, "%3s%3s", "#", "");
  fprintf (fp, "%20s", "x");
  for (i = 0; i < SIZE; ++i)
    {
      fprintf (fp, "%20d", i);
    }
  fprintf (fp, "\n");

  fprintf (fp, "%3s%3s", "#", "");
  fprintf (fp, "%20s", "--------------------");
  for (i = 0; i < SIZE; ++i)
    {
      fprintf (fp, "%20s", "--------------------");
    }
  fprintf (fp, "\n");

  for (i = 0; i < SIZE_HIGH_RES; ++i)
    {
      fprintf (fp, "%6s", "");
      fprintf (fp, "%20.4e", high_res_ordinate[i]);
      for (j = 0; j < SIZE; ++j)
        {
          fprintf (fp, "%20.4e",
                   rhs[j] * chebyshev_poly_1 (j, high_res_ordinate[i]));
        }
      fprintf (fp, "\n");
    }

  fclose (fp);
  printf ("Saved basis function data to file '%s'.\n",
          basis_functions_filename);

  /* Write the basis functions to a file. */
  /* TODO: make the file name an adjustable parameter. */

  fp = fopen (coeffs_filename, "w");
  if (fp == NULL)
    {
      fprintf (stderr, "Can't open file '%s' for writing!\n",
               coeffs_filename);
      return -1;
    }

  for (i = 0; i < SIZE; ++i)
    {
      fprintf (fp, "%10d %20.4e\n", i, rhs[i]);
    }

  fclose (fp);
  printf ("Saved coefficient data to file '%s'.\n", coeffs_filename);


  return 0;
}
