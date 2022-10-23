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
    Vector biased() const;
    explicit Vector(ifstream& inputFile);
    void add(double x);
    void insert(int pos, double x);
    void remove(int pos);
    void clear();
    double sumOfElements() const;
    unsigned long maxIndex() const;
    void sigmoid();
    void tanh();
    void relu();
    void reluDerivative();
    Vector skipVector(unsigned long mod, unsigned long value) const;
    void add(const Vector& v);
    void subtract(const Vector& v);
    Vector difference(const Vector& v) const;
    double dotProduct(const Vector& v) const;
    double dotProduct() const;
    Vector elementProduct(const Vector& v) const;
    void divide(double value);
    void multiply(double value);
    Vector product(double value) const;
    void l1Normalize();
    double l2Norm() const;
    unsigned long getSize() const;
    double cosineSimilarity(const Vector& v) const;
    double getValue(unsigned long index) const;
    void setValue(unsigned long index, double value);
    void addValue(unsigned long index, double value);
    double sum() const;
    void swap(int index1, int index2);
    void serialize(ostream& outputFile);
};


#endif //MATH_VECTOR_H
