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
    Matrix(int row, int col);
    Matrix(int row, int col, double min, double max, default_random_engine randomEngine);
    explicit Matrix(int size);
    Matrix(const Vector& v1, const Vector& v2);
    Matrix clone();
    void printToFile(const string& fileName) const;
    double getValue(int rowNo, int colNo) const;
    void setValue(int rowNo, int colNo, double value);
    void addValue(int rowNo, int colNo, double value);
    void increment(int rowNo, int colNo);
    int getRow() const;
    Vector getRow(int row) const;
    vector<double> getColumn(int column) const;
    int getColumn() const;
    void columnWiseNormalize();
    void multiplyWithConstant(double constant);
    void divideByConstant(double constant);
    void add(const Matrix& m);
    void add(int rowNo, const Vector& v);
    void subtract(const Matrix& m);
    Vector multiplyWithVectorFromLeft(const Vector& v) const;
    Vector multiplyWithVectorFromRight(const Vector& v) const;
    double columnSum(int columnNo) const;
    Vector sumOfRows() const;
    double rowSum(int rowNo) const;
    Matrix multiply(const Matrix& m) const;
    Matrix elementProduct(const Matrix& m) const;
    double sumOfElements() const;
    double trace() const;
    Matrix transpose();
    Matrix partial(int rowstart, int rowend, int colstart, int colend);
    bool isSymmetric();
    double determinant();
    void inverse();
    Matrix choleskyDecomposition();
    vector<Eigenvector> characteristics();
private:
    void rotate(double s, double tau, int i, int j, int k, int l);
};


#endif //MATH_MATRIX_H
