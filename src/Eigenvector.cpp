#include <utility>

//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#include "Eigenvector.h"

/**
 * A constructor of Eigenvector which takes a double eigenValue and an vector values as inputs.
 * It calls its super class Vector with values vector and initializes eigenValue variable with its
 * eigenValue input.
 *
 * @param eigenValue double input.
 * @param values     vector input.
 */
Eigenvector::Eigenvector(double eigenValue, const vector<double>& values) : Vector(values) {
    this->eigenValue = eigenValue;
}

/**
 * The getEigenValue method which returns the eigenValue variable.
 *
 * @return eigenValue variable.
 */
double Eigenvector::getEigenValue() const{
    return eigenValue;
}
