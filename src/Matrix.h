//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H
#include<string>
#include<random>
#include "Vector.h"
#include "Eigenvector.h"

using namespace std;

class Matrix {
private:
    int row;
    int col;
    double** values;
public:
    explicit Matrix(const string& filename);
    explicit Matrix(ifstream& inputStream);
    Matrix(int row, int col);
    Matrix(int row, int col, double min, double max, default_random_engine randomEngine);
    explicit Matrix(int size);
    Matrix(const Vector& v1, const Vector& v2);
    [[nodiscard]] Matrix clone() const;
    ~Matrix();
    void printToFile(const string& fileName) const;
    [[nodiscard]] double getValue(int rowNo, int colNo) const;
    void setValue(int rowNo, int colNo, double value) const;
    void addValue(int rowNo, int colNo, double value) const;
    void increment(int rowNo, int colNo) const;
    [[nodiscard]] int getRow() const;
    [[nodiscard]] Vector getRow(int row) const;
    [[nodiscard]] vector<double> getColumn(int column) const;
    [[nodiscard]] int getColumn() const;
    void columnWiseNormalize() const;
    void multiplyWithConstant(double constant) const;
    void divideByConstant(double constant) const;
    void add(const Matrix& m) const;
    [[nodiscard]] Matrix sum(const Matrix& m) const;
    void add(int rowNo, const Vector& v) const;
    void subtract(const Matrix& m) const;
    [[nodiscard]] Matrix difference(const Matrix& m) const;
    [[nodiscard]] Vector multiplyWithVectorFromLeft(const Vector& v) const;
    [[nodiscard]] Vector multiplyWithVectorFromRight(const Vector& v) const;
    [[nodiscard]] double columnSum(int columnNo) const;
    [[nodiscard]] Vector sumOfRows() const;
    [[nodiscard]] double rowSum(int rowNo) const;
    [[nodiscard]] Matrix multiply(const Matrix& m) const;
    [[nodiscard]] Matrix elementProduct(const Matrix& m) const;
    [[nodiscard]] double sumOfElements() const;
    [[nodiscard]] double trace() const;
    [[nodiscard]] Matrix transpose() const;
    [[nodiscard]] Matrix partial(int rowstart, int rowend, int colstart, int colend) const;
    [[nodiscard]] bool isSymmetric() const;
    [[nodiscard]] double determinant() const;
    void inverse() const;
    [[nodiscard]] Matrix choleskyDecomposition() const;
    [[nodiscard]] vector<Eigenvector> characteristics() const;
private:
    void rotate(double s, double tau, int i, int j, int k, int l) const;
};


#endif //MATH_MATRIX_H
