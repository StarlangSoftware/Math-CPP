cmake_minimum_required(VERSION 3.12)
project(Math)

set(CMAKE_CXX_STANDARD 20)

add_library(Math Distribution.cpp Distribution.h DeterminantZero.cpp DeterminantZero.h VectorSizeMismatch.cpp VectorSizeMismatch.h MatrixColumnMismatch.cpp MatrixColumnMismatch.h MatrixDimensionMismatch.cpp MatrixDimensionMismatch.h MatrixNotPositiveDefinite.cpp MatrixNotPositiveDefinite.h MatrixNotSymmetric.cpp MatrixNotSymmetric.h MatrixRowColumnMismatch.cpp MatrixRowColumnMismatch.h MatrixRowMismatch.cpp MatrixRowMismatch.h Vector.cpp Vector.h Eigenvector.cpp Eigenvector.h DiscreteDistribution.cpp DiscreteDistribution.h Matrix.cpp Matrix.h MatrixNotSquare.cpp MatrixNotSquare.h)

add_executable(TestMath Distribution.cpp Distribution.h DeterminantZero.cpp DeterminantZero.h VectorSizeMismatch.cpp VectorSizeMismatch.h MatrixColumnMismatch.cpp MatrixColumnMismatch.h MatrixDimensionMismatch.cpp MatrixDimensionMismatch.h MatrixNotPositiveDefinite.cpp MatrixNotPositiveDefinite.h MatrixNotSymmetric.cpp MatrixNotSymmetric.h MatrixRowColumnMismatch.cpp MatrixRowColumnMismatch.h MatrixRowMismatch.cpp MatrixRowMismatch.h Vector.cpp Vector.h Eigenvector.cpp Eigenvector.h DiscreteDistribution.cpp DiscreteDistribution.h Matrix.cpp Matrix.h MatrixNotSquare.cpp MatrixNotSquare.h TestMath.cpp)