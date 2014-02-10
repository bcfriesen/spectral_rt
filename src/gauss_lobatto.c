#include <math.h>
#include <assert.h>

#include <globals.h>

double gauss_lobatto (const int i) {
    double pi;

    assert (i >= 0 && i < n);
    pi = acos(-1.0);
    return cos(pi*(double)i/(double)(n-1));
}
