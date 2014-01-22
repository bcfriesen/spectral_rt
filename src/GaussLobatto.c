#include <math.h>

// Returns the kth Gauss-Lobatto interpolation point
double gauss_lobatto_point(int k, int N)
{
    double pi = acos(-1.0);
    
    double x_k = cos(k * pi / N);

    return x_k;
}
