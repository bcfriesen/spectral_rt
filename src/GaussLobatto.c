#include <math.h>
#include <stdio.h>

// Returns the kth Gauss-Lobatto interpolation point
double gauss_lobatto_point(int k, int N)
{
    double x_k = - cos(k * M_PI / N);

    return x_k;
}

int main(int argc, char *argv[])
{
    int k = 1;
    int N = 4;

    double test_point = gauss_lobatto_point(k, N);

    printf("test_point = %f\n", test_point);

    return 0;
}
