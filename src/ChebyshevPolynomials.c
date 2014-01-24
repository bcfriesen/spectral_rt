#include <math.h>

#include <cctk.h>

double cheby_poly_1(int n, double x)
{
    int i;
    double T_n, T_nm1, T_nm2;

    if (n < 0) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Invalid Chebyshev polynomial (1st kind) index: %d", n);
    } else if (n == 0) {
        return (1.0);
    } else if (n == 1) {
        return (x);
    } else {
        T_nm2 = 1.0;
        T_nm1 = x;
        for (i = 2; i <= n; i++) {
            T_n = 2.0*x*T_nm1 - T_nm2;
            T_nm2 = T_nm1;
            T_nm1 = T_n;
        }
        return (T_n);
    }
}

double cheby_poly_2(int n, double x)
{
    int i;
    double U_n, U_nm1, U_nm2;
    if (n < 0) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Invalid Chebyshev polynomial (2nd kind) index: %d", n);
    } else if (n == 0) {
        return (1.0);
    } else if (n == 1) {
        return (2.0 * x);
    } else {
        U_nm2 = 1.0;
        U_nm1 = 2.0 * x;
        for (i = 2; i <= n; i++) {
            U_n = 2.0 * x * U_nm1 - U_nm2;
            U_nm2 = U_nm1;
            U_nm1 = U_n;
        }
        return (U_n);
    }
}
