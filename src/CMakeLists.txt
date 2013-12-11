set (CMAKE_C_FLAGS "-Wall")
set (CMAKE_EXE_LINKER_FLAGS "-Wall")

include_directories("/usr/local/include"
                    ".")

set (test_galerkin_method_SOURCES
     approx_result.c
     chebyshev_polynomials.c
     main.c
    )

ADD_EXECUTABLE(test_galerkin_method ${test_galerkin_method_SOURCES})

target_link_libraries(test_galerkin_method m gsl gslcblas lapack blas)
