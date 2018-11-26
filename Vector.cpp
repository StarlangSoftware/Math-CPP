//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "Vector.h"

/**
 * A constructor of {@link Vector} class which takes an {@link vector} values as an input. Then, initializes
 * values {@link vector} and size variable with given input and ts size.
 *
 * @param values {@link vector} input.
 */
Vector::Vector(vector<double> values) {
    this->values = values;
    size = values.capacity();
}

/**
 * Another constructor of {@link Vector} class which takes integer size and double x as inputs. Then, initializes size
 * variable with given size input and creates new values {@link vector} and adds given input x to values {@link vector}.
 *
 * @param size {@link vector} size.
 * @param x    item to add values {@link vector}.
 */
Vector::Vector(unsigned long size, double x) {
    this->size = size;
    values.assign(size, x);
}

/**
 * Another constructor of {@link Vector} class which takes integer size, integer index and double x as inputs. Then, initializes size
 * variable with given size input and creates new values {@link vector} and adds 0.0 to values {@link vector}.
 * Then, sets the item of values {@link vector} at given index as given input x.
 *
 * @param size  {@link vector} size.
 * @param index to set a particular item.
 * @param x     item to add values {@link vector}'s given index.
 */
Vector::Vector(unsigned long size, int index, double x) {
    this->size = size;
    values.assign(size, 0.0);
    values[index] = x;
}
