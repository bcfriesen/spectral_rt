cmake_minimum_required (VERSION 2.8)
project (test_galerkin_method)
enable_language (C)

set (CMAKE_C_COMPILER "clang")

set (CMAKE_VERBOSE_MAKEFILE "true")

add_subdirectory (src)
