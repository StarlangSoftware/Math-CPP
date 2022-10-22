//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H
#include <vector>
using namespace std;

class Vector {
private:
    unsigned long size;
    vector<double> values;
public:
    explicit Vector(const vector<double>& values);
    Vector(unsigned long size, double x);
    Vector(unsigned long size, int index, double x);
    Vector(double* values, unsigned long size);
    Vector biased();
    explicit Vector(ifstream& inputFile);
    void add(double x);
    void insert(int pos, double x);
    void remove(int pos);
    void clear();
    double sumOfElements();
    unsigned long maxIndex();
    void sigmoid();
    void tanh();
    void relu();
    void reluDerivative();
    Vector skipVector(unsigned long mod, unsigned long value);
    void add(const Vector& v);
    void subtract(const Vector& v);
    Vector difference(const Vector& v);
    double dotProduct(const Vector& v);
    double dotProduct();
    Vector elementProduct(const Vector& v);
    void divide(double value);
    void multiply(double value);
    Vector product(double value);
    void l1Normalize();
    double l2Norm() const;
    unsigned long getSize() const;
    double cosineSimilarity(const Vector& v);
    double getValue(unsigned long index) const;
    void setValue(unsigned long index, double value);
    void addValue(unsigned long index, double value);
    double sum();
    void swap(int index1, int index2);
    void serialize(ostream& outputFile);
};


#endif //MATH_VECTOR_H
