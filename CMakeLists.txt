cmake_minimum_required(VERSION 3.12)
project(Math)

set(CMAKE_CXX_STANDARD 17)

add_library(Math Distribution.cpp Distribution.h DeterminantZero.cpp DeterminantZero.h VectorSizeMismatch.cpp VectorSizeMismatch.h MatrixColumnMismatch.cpp MatrixColumnMismatch.h MatrixDimensionMismatch.cpp MatrixDimensionMismatch.h MatrixNotPositiveDefinite.cpp MatrixNotPositiveDefinite.h MatrixNotSymmetric.cpp MatrixNotSymmetric.h MatrixRowColumnMismatch.cpp MatrixRowColumnMismatch.h MatrixRowMismatch.cpp MatrixRowMismatch.h Vector.cpp Vector.h)