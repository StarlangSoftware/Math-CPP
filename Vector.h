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
    Vector(vector<double> values);
    Vector(unsigned long size, double x);
    Vector(unsigned long size, int index, double x);
};


#endif //MATH_VECTOR_H
