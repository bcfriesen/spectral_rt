#include  <error_codes.h>

double
chebyshev_poly_1 (const int n, const double x, int *err)
{

  double T_n, T_nm1, T_nm2;
  int i;

  if (n == 0)
    {
      T_n = 1.0;
    }
  else if (n == 1)
    {
      T_n = x;
    }
  else if (n > 1)
    {
      T_nm1 = x;
      T_nm2 = 1.0;
      for (i = 2; i <= n; ++i)
	{
	  T_n = 2.0 * x * T_nm1 - T_nm2;
	  T_nm2 = T_nm1;
	  T_nm1 = T_n;
	}
    }
  else
    {
      T_n = 0.0;
      *err = ENEGNCHEBY;
    }

  if (n >= 0)
    *err = SUCCESS;

  return (T_n);
}

double
chebyshev_poly_2 (const int n, const double x, int *err)
{

  double U_n, U_nm1, U_nm2;
  int i;

  if (n == 0)
    {
      U_n = 1.0;
    }
  else if (n == 1)
    {
      U_n = 2.0 * x;
    }
  else if (n > 1)
    {
      U_nm1 = 2.0 * x;
      U_nm2 = 1.0;
      for (i = 2; i <= n; ++i)
	{
	  U_n = 2.0 * x * U_nm1 - U_nm2;
	  U_nm2 = U_nm1;
	  U_nm1 = U_n;
	}
    }
  else
    {
      U_n = 0.0;
      *err = ENEGNCHEBY;
    }

  if (n >= 0)
    *err = SUCCESS;

  return (U_n);
}
