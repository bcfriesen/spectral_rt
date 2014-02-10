#include <assert.h>

#include <globals.h>

int c_bar(const int j) {
    int return_val;
    assert(j >= 0 && j < n);

    if (j == 0 || j == (n-1)) {
        return_val = 2;
    } else {
        return_val = 1;
    }

    return (return_val);
}
