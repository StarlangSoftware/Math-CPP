#include <utility>

//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#include "Eigenvector.h"

/**
 * A constructor of {@link Eigenvector} which takes a double eigenValue and an {@link vector} values as inputs.
 * It calls its super class {@link Vector} with values {@link vector} and initializes eigenValue variable with its
 * eigenValue input.
 *
 * @param eigenValue double input.
 * @param values     {@link vector} input.
 */
Eigenvector::Eigenvector(double eigenValue, vector<double> values) : Vector(move(values)) {
    this->eigenValue = eigenValue;
}

/**
 * The getEigenValue method which returns the eigenValue variable.
 *
 * @return eigenValue variable.
 */
double Eigenvector::getEigenValue() {
    return eigenValue;
}
