#include <math.h>

#include <globals.h>
#include <c_bar.h>

void calc_chebyshev_derivative_matrix(const double *grid, double **chebyshev_derivative_matrix) {
    int j, l;
    for (l=0; l<n; l++) {
        for (j=0; j<n; j++) {
            if (j != l) {
                chebyshev_derivative_matrix[j][l] = ((double)c_bar(j) * pow(-1.0,j+l)) / ((double)c_bar(l) * (grid[j]-grid[l]));
            } else {
                if (j >= 1 && j <= n-2) {
                    chebyshev_derivative_matrix[j][l] = -grid[l] / (2.0 * (1.0 - pow(grid[l],2)));
                } else if (j == 0) {
                    chebyshev_derivative_matrix[j][l] = (2.0 * pow((double)(n-1),2) + 1.0) / 6.0;
                } else if (j == n-1) {
                    chebyshev_derivative_matrix[j][l] = -(2.0 * pow((double)(n-1),2) + 1.0) / 6.0;
                }
            }
        }
    }
}
