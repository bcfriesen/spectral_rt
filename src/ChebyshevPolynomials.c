#include <math.h>

double cheby_poly_1(unsigned int n, double x)
{
    return (cos((double)n * acos(x)));
}

double cheby_poly_2(unsigned int n, double x)
{
    return (sin((double)(n+1) * acos(x)) / sin(acos(x)));
}
