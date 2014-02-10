include_directories(".")

set (test_spectral_rt_SOURCES
     c_bar.c
     chebyshev_polynomials.c
     gauss_lobatto.c
     main.c
    )

ADD_EXECUTABLE(test_spectral_rt ${test_spectral_rt_SOURCES})

target_link_libraries(test_spectral_rt m lapack blas fftw3)
