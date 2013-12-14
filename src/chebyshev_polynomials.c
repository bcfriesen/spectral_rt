#include <assert.h>
#include <math.h>

#include  <error_codes.h>

/* TODO: figure out at what n it becomes faster to use the trigonometric
 * definitions rather than the recurrence relations. */

double
chebyshev_poly_1 (const int n, const double x)
{
  assert (n >= 0);
  return (cos ((double) n * acos (x)));
}

double
chebyshev_poly_2 (const int n, const double x)
{
  assert (n >= 0);
  return (sin ((n + 1) * acos (x)) / sin (acos (x)));
}
