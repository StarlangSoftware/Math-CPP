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
    vector<Vector> values;
public:
    explicit Matrix(const string& filename);
    Matrix(int row, int col);
    Matrix(int row, int col, double min, double max, default_random_engine randomEngine);
    explicit Matrix(int size);
    Matrix(const Vector& v1, const Vector& v2);
    explicit Matrix(ifstream& inputFile);
    Matrix clone();
    void printToFile(const string& fileName);
    double getValue(int rowNo, int colNo);
    void setValue(int rowNo, int colNo, double value);
    void addValue(int rowNo, int colNo, double value);
    void increment(int rowNo, int colNo);
    int getRow();
    Vector getRow(int row);
    vector<double> getColumn(int column);
    int getColumn();
    void columnWiseNormalize();
    void multiplyWithConstant(double constant);
    void divideByConstant(double constant);
    void add(const Matrix& m);
    void add(int rowNo, const Vector& v);
    void subtract(const Matrix& m);
    Vector multiplyWithVectorFromLeft(const Vector& v);
    Vector multiplyWithVectorFromRight(const Vector& v);
    double columnSum(int columnNo);
    Vector sumOfRows();
    double rowSum(int rowNo);
    Matrix multiply(const Matrix& m);
    Matrix elementProduct(const Matrix& m);
    double sumOfElements();
    double trace();
    Matrix transpose();
    Matrix partial(int rowstart, int rowend, int colstart, int colend);
    bool isSymmetric();
    double determinant();
    void inverse();
    Matrix choleskyDecomposition();
    vector<Eigenvector> characteristics();
    void serialize(ostream& outputFile);
private:
    void rotate(double s, double tau, int i, int j, int k, int l);
};


#endif //MATH_MATRIX_H
