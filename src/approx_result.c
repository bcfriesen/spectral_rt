#include <chebyshev_polynomials.h>

double
approx_result (const double x, const double coeffs[], const int SIZE)
{

  int i;
  double result;

  result = 0.0;

  for (i = 0; i < SIZE; ++i)
    {
      result = result + coeffs[i] * chebyshev_poly_1 (i, x);
    }

  return result;
}
