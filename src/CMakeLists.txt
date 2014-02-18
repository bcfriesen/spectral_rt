include_directories("."
    "/home/brian/src/lapack-3.5.0/build/include")

find_library(LAPACKE_LIBS lapacke)
find_library(LAPACK_LIBS lapack)
find_library(BLAS_LIBS blas)

set (test_spectral_rt_SOURCES
     c_bar.c
     calc_chebyshev_derivative_matrix_1.c
     chebyshev_polynomials.c
     gauss_lobatto.c
     main.c
    )

ADD_EXECUTABLE(test_spectral_rt ${test_spectral_rt_SOURCES})

target_link_libraries(test_spectral_rt m fftw3 ${LAPACKE_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} gfortran)
