//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#ifndef MATH_EIGENVECTOR_H
#define MATH_EIGENVECTOR_H
#include "Vector.h"

class Eigenvector : public Vector{
private:
    double eigenValue;
public:
    Eigenvector(double eigenValue, const vector<double>& values);
    [[nodiscard]] double getEigenValue() const;
};


#endif //MATH_EIGENVECTOR_H
