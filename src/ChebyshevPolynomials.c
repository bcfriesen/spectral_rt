#include <math.h>

#include <cctk.h>

double cheby_poly_1(int n, double x)
{
    if (n < 0){
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Invalid Chebyshev polynomial (1st kind) index: %d", n);
    }
    return (cos((double)n * acos(x)));
}

double cheby_poly_2(int n, double x)
{
    if (n < 0){
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Invalid Chebyshev polynomial (2nd kind) index: %d", n);
    }
    return (sin((double)(n+1) * acos(x)) / sin(acos(x)));
}
