//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H
#include <vector>
#include <fstream>
using namespace std;

class Vector {
private:
    unsigned long size;
    vector<double> values;
public:
    explicit Vector(const vector<double>& values);
    Vector(unsigned long size, double x);
    Vector(unsigned long size, int index, double x);
    Vector(const double* values, unsigned long size);
    [[nodiscard]] Vector biased() const;
    explicit Vector(ifstream& inputFile);
    void add(double x);
    void insert(int pos, double x);
    void remove(int pos);
    void clear();
    [[nodiscard]] double sumOfElements() const;
    [[nodiscard]] unsigned long maxIndex() const;
    void sigmoid();
    void tanh();
    void relu();
    void reluDerivative();
    [[nodiscard]] Vector skipVector(unsigned long mod, unsigned long value) const;
    void add(const Vector& v);
    void subtract(const Vector& v);
    [[nodiscard]] Vector difference(const Vector& v) const;
    [[nodiscard]] double dotProduct(const Vector& v) const;
    [[nodiscard]] double dotProduct() const;
    [[nodiscard]] Vector elementProduct(const Vector& v) const;
    void divide(double value);
    void multiply(double value);
    [[nodiscard]] Vector product(double value) const;
    void l1Normalize();
    [[nodiscard]] double l2Norm() const;
    [[nodiscard]] unsigned long getSize() const;
    [[nodiscard]] double cosineSimilarity(const Vector& v) const;
    [[nodiscard]] double getValue(unsigned long index) const;
    void setValue(unsigned long index, double value);
    void addValue(unsigned long index, double value);
    [[nodiscard]] double sum() const;
    void swap(int index1, int index2);
    void serialize(ostream& outputFile) const;
};


#endif //MATH_VECTOR_H
