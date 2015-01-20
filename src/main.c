#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <fftw3.h>

#include <chebyshev_polynomials.h>
#include <globals.h>

const int n = 128;
double *grid;

int main () {
  int i;
  double *spec;
  double *spec_deriv;
  double *x_prime;
  double *f_hat_prime;
  double *f_hat;
  double *f_tilde_prime;
  double *f_tilde;
  double *F_tilde;
  double *F;
  double *f_n;
  fftw_plan plan;
  double tmp;
  FILE *fp;
  double *c_bar;
  double pi;

  pi = acos(-1.0);

  printf("malloc-ing arrays ... ");
  spec = malloc(n * sizeof(double));
  spec_deriv = malloc(n * sizeof(double));
  grid = malloc(n * sizeof(double));
  x_prime = malloc(n * sizeof(double));
  f_hat_prime = malloc(n * sizeof(double));
  f_hat = malloc(n * sizeof(double));
  f_tilde_prime = malloc(n * sizeof(double));
  f_tilde = malloc(n * sizeof(double));
  F_tilde = malloc(n * sizeof(double));
  F = malloc(n * sizeof(double));
  f_n = malloc(n * sizeof(double));
  c_bar = malloc(n * sizeof(double));
  printf("done!\n");


  printf("setting up Gauss-Lobatto grid ... ");
  /* Set up ordinate grid using Gauss-Lobatto points. */
  for (i = 0; i < n; i++) {
      grid[i] = cos(pi*(double)i/(double)(n-1));
  }
  printf("done!\n");

  /* Fill in c_bar array. */
  c_bar[0] = 2.0;
  c_bar[n-1] = 2.0;
  for (i = 1; i < n-1; i++) {
    c_bar[i] = 1.0;
  }

  for (i = 0; i < n; i++) {
      x_prime[i] = grid[i] / c_bar[i];
  }

  /* This FFT gives us f_hat_prime. */
  printf("Setting up FFT to get f_hat_prime ... ");
  plan = fftw_plan_r2r_1d(n, x_prime, f_hat_prime, FFTW_REDFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  printf("done!\n");

  /* From f_hat_prime calculate f_hat. */
  for (i = 0; i < n; i++) {
      f_hat[i] = f_hat_prime[i] + x_prime[0] + pow(-1.0, i)*x_prime[n-1];
  }

  /* From f_hat calculate f_tilde_prime. */
  for (i = 0; i < n; i++) {
      f_tilde_prime[i] = (2.0 * f_hat[i]) / ((double)(n-1) * c_bar[i]);
  }

  /* From f_tilde_prime calcualte all f_tilde except f_tilde[0]. */
  f_tilde[n-1] = (1.0 / (double)(2*(n-1))) * c_bar[n-2]*f_tilde_prime[n-2]; /* Do this one separately because it goes "out of bounds" of f_tilde_prime. */
  for (i = n-2; i >= 1; i--) {
      f_tilde[i] = (1.0 / (double)(2*i)) * (c_bar[i-1]*f_tilde_prime[i-1] - f_tilde_prime[i+1]);
  }
  /* Calculate f_tilde[0] by using the boundary condition. */
  tmp = 0.0;
  for (i = 1; i < n; i++) {
      tmp += f_tilde[i] * chebyshev_poly_1(i, 0.0);
  }
  f_tilde[0] = -tmp;

  /* From f_tilde calculate F_tilde. */
  for (i = 0; i < n; i++) {
      F_tilde[i] = f_tilde[i] / 2.0;
  }

  /* This FFT gives us F. */
  printf("Setting up FFT to get F ... ");
  plan = fftw_plan_r2r_1d(n, F_tilde, F, FFTW_REDFT00, FFTW_ESTIMATE);
  fftw_execute(plan);
  printf("done!\n");

  printf("Writing results to file ... ");
  /* From F calculate f_n. */
  fp = fopen("results.dat", "w");
  fprintf(fp, "%1s %15s %15s %15s\n", "#", "x[i]", "exact", "approx");
  for (i = 0; i < n; i++) {
      f_n[i] = F[i] + F_tilde[0] + pow(-1.0,i)*F_tilde[n-1];
      fprintf(fp, "%1s %+15.7e %+15.7e %+15.7e\n", "", grid[i], pow(grid[i],2), f_n[i]);
  }
  printf("done!\n");

  free(grid);
  free(spec);
  free(spec_deriv);
  free(x_prime);
  free(f_hat_prime);
  free(f_hat);
  free(f_tilde_prime);
  free(f_tilde);
  free(F_tilde);
  free(F);
  free(f_n);
  fftw_free(plan);
  fclose(fp);
  free(c_bar);

  return 0;
}
