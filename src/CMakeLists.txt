include_directories(".")

set (test_spectral_rt_SOURCES
     chebyshev_polynomials.c
     main.c
    )

ADD_EXECUTABLE(test_spectral_rt ${test_spectral_rt_SOURCES})

target_link_libraries(test_spectral_rt m fftw3)
