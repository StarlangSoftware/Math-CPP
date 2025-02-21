//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include <cmath>
#include <fstream>
#include "Vector.h"
#include "VectorSizeMismatch.h"

/**
 * A constructor of Vector class which takes an vector values as an input. Then, initializes
 * values vector and size variable with given input and ts size.
 *
 * @param values vector input.
 */
Vector::Vector(const vector<double>& values) {
    this->values = values;
    size = values.size();
}

/**
 * Another constructor of Vector class which takes integer size and double x as inputs. Then, initializes size
 * variable with given size input and creates new values vector and adds given input x to values vector.
 *
 * @param size vector size.
 * @param x    item to add values vector.
 */
Vector::Vector(unsigned long size, double x) {
    this->size = size;
    values.assign(size, x);
}

/**
 * Another constructor of Vector class which takes integer size, integer index and double x as inputs. Then, initializes size
 * variable with given size input and creates new values vector and adds 0.0 to values vector.
 * Then, sets the item of values vector at given index as given input x.
 *
 * @param size  vector size.
 * @param index to set a particular item.
 * @param x     item to add values vector's given index.
 */
Vector::Vector(unsigned long size, int index, double x) {
    this->size = size;
    values.assign(size, 0.0);
    values[index] = x;
}

/**
 * Another constructor of Vector class which takes double values array as an input.
 * It creates new values vector and adds given input values array's each item to the values vector.
 * Then, initializes size with given values input array's length.
 *
 * @param values double array input.
 */
Vector::Vector(double* values, unsigned long size) {
    this->values.reserve(size);
    for (int i = 0; i < size; i++) {
        this->values.push_back(values[i]);
    }
    this->size = size;
}

/**
 * The biased method creates a Vector result, add adds each item of values vector into the result Vector.
 * Then, insert 1.0 to 0th position and return result Vector.
 *
 * @return result Vector.
 */
Vector Vector::biased() const{
    Vector result = Vector(0, 0.0);
    for (double value : values) {
        result.add(value);
    }
    result.insert(0, 1.0);
    return result;
}

/**
 * The add method adds given input to the values vector and increments the size variable by one.
 *
 * @param x double input to add values vector.
 */
void Vector::add(double x) {
    values.push_back(x);
    size++;
}

/**
 * The insert method puts given input to the given index of values vector and increments the size variable by one.
 *
 * @param pos index to insert input.
 * @param x   input to insert to given index of values vector.
 */
void Vector::insert(int pos, double x) {
    values.insert(values.begin() + pos, x);
    size++;
}

/**
 * The remove method deletes the item at given input position of values vector and decrements the size variable by one.
 *
 * @param pos index to remove from values vector.
 */
void Vector::remove(int pos) {
    values.erase(values.begin() + pos);
    size--;
}

/**
 * The clear method sets all the elements of values vector to 0.0.
 */
void Vector::clear() {
    for (double & value : values) {
        value = 0.0;
    }
}

/**
 * The sumOfElements method sums up all elements in the vector.
 *
 * @return Sum of all elements in the vector.
 */
double Vector::sumOfElements() const{
    double total = 0;
    for (int i = 0; i < size; i++) {
        total += values[i];
    }
    return total;
}

/**
 * The maxIndex method gets the first item of values ArrayList as maximum item, then it loops through the indices
 * and if a greater value than the current maximum item comes, it updates the maximum item and returns the final
 * maximum item's index.
 *
 * @return final maximum item's index.
 */
unsigned long Vector::maxIndex() const{
    unsigned long index = 0;
    double max = values.at(0);
    for (unsigned long i = 1; i < size; i++) {
        if (values.at(i) > max) {
            max = values.at(i);
            index = i;
        }
    }
    return index;
}

/**
 * The sigmoid method loops through the values vector and sets each ith item with sigmoid function, i.e
 * 1 / (1 + exp(-values.get(i))), i ranges from 0 to size.
 */
void Vector::sigmoid() {
    for (unsigned long i = 0; i < size; i++) {
        values[i] = 1 / (1 + exp(-values.at(i)));
    }
}

/**
 * The tanh method loops through the values ArrayList and sets each ith item with tanh function.
 */
void Vector::tanh() {
    for (unsigned long i = 0; i < size; i++) {
        values[i] = ::tanh(values[i]);
    }
}

/**
 * The relu method loops through the values ArrayList and sets each ith item with relu function.
 */
void Vector::relu() {
    for (unsigned long i = 0; i < size; i++) {
        if (values[i] < 0){
            values[i] = 0.0;
        }
    }
}

/**
 * The reluDerivative method loops through the values ArrayList and sets each ith item with the derivative of
 * relu function.
 */
void Vector::reluDerivative() {
    for (unsigned long i = 0; i < size; i++) {
        if (values[i] > 0){
            values[i] = 1.0;
        } else {
            values[i] = 0.0;
        }
    }
}

/**
 * The skipVector method takes a mod and a value as inputs. It creates a new result Vector, and assigns given input value to i.
 * While i is less than the size, it adds the ith item of values vector to the result and increments i by given mod input.
 *
 * @param mod   integer input.
 * @param value integer input.
 * @return result Vector.
 */
Vector Vector::skipVector(unsigned long mod, unsigned long value) const{
    Vector result = Vector(0, 0.0);
    unsigned long i = value;
    while (i < size) {
        result.add(values.at(i));
        i += mod;
    }
    return result;
}

/**
 * The add method takes a Vector v as an input. It sums up the corresponding elements of both given vector's
 * values vector and values vector and puts result back to the values vector.
 * If their sizes do not match, it throws a VectorSizeMismatch exception.
 *
 * @param v Vector to add.
 */
void Vector::add(const Vector& v) {
    if (size != v.size) {
        throw VectorSizeMismatch();
    }
    for (int i = 0; i < size; i++) {
        values[i] = values.at(i) + v.values.at(i);
    }
}

/**
 * The subtract method takes a Vector v as an input. It subtracts the corresponding elements of given vector's
 * values vector from values vector and puts result back to the values vector.
 * If their sizes do not match, it throws a VectorSizeMismatch exception.
 *
 * @param v Vector to subtract from values vector.
 */
void Vector::subtract(const Vector& v) {
    if (size != v.size) {
        throw VectorSizeMismatch();
    }
    for (unsigned long i = 0; i < size; i++) {
        values[i] = values.at(i) - v.values.at(i);
    }
}

/**
 * The difference method takes a Vector v as an input. It creates a new double array result, then
 * subtracts the corresponding elements of given vector's values vector from values vector and puts
 * result back to the result array. If their sizes do not match, it throws a VectorSizeMismatch exception.
 *
 * @param v Vector to find difference from values vector.
 * @return new Vector with result array.
 */
Vector Vector::difference(const Vector& v) const{
    if (size != v.size) {
        throw VectorSizeMismatch();
    }
    vector<double> result(size);
    for (unsigned long i = 0; i < size; i++) {
        result[i] = values.at(i) - v.values.at(i);
    }
    return Vector(result);
}

/**
 * The dotProduct method takes a Vector v as an input. It creates a new double variable result, then
 * multiplies the corresponding elements of given vector's values vector with values vector and assigns
 * the multiplication to the result. If their sizes do not match, it throws a VectorSizeMismatch exception.
 *
 * @param v Vector to find dot product.
 * @return double result.
 */
double Vector::dotProduct(const Vector& v) const{
    if (size != v.size) {
        throw VectorSizeMismatch();
    }
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += values.at(i) * v.values.at(i);
    }
    return result;
}

/**
 * The dotProduct method creates a new double variable result, then squares the elements of values vector and assigns
 * the accumulation to the result.
 *
 * @return double result.
 */
double Vector::dotProduct() const{
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += values.at(i) * values.at(i);
    }
    return result;
}

/**
 * The elementProduct method takes a Vector v as an input. It creates a new double array result, then
 * multiplies the corresponding elements of given vector's values vector with values ArrayList and assigns
 * the multiplication to the result array. If their sizes do not match, it throws a VectorSizeMismatch exception.
 *
 * @param v Vector to find dot product.
 * @return Vector with result array.
 */
Vector Vector::elementProduct(const Vector& v) const{
    if (size != v.size) {
        throw VectorSizeMismatch();
    }
    vector<double> result(size);
    for (unsigned long i = 0; i < size; i++) {
        result[i] = values.at(i) * v.values.at(i);
    }
    return Vector(result);
}

/**
 * The divide method takes a double value as an input and divides each item of values vector with given value.
 *
 * @param value is used to divide items of values vector.
 */
void Vector::divide(double value) {
    for (unsigned long i = 0; i < size; i++) {
        values[i] = values.at(i) / value;
    }
}

/**
 * The multiply method takes a double value as an input and multiplies each item of values vector with given value.
 *
 * @param value is used to multiply items of values vector.
 */
void Vector::multiply(double value) {
    for (unsigned long i = 0; i < size; i++) {
        values[i] = values.at(i) * value;
    }
}

/**
 * The product method takes a double value as an input and creates a new result Vector, then multiplies each
 * item of values vector with given value and adds to the result Vector.
 *
 * @param value is used to multiply items of values vector.
 * @return Vector result.
 */
Vector Vector::product(double value) const{
    Vector result = Vector(0, 0.0);
    for (unsigned long i = 0; i < size; i++) {
        result.add(values.at(i) * value);
    }
    return result;
}

/**
 * The l1Normalize method is used to apply Least Absolute Errors, it accumulates items of values vector and sets
 * each item by dividing it by the summation value.
 */
void Vector::l1Normalize() {
    double sum = 0;
    for (unsigned long i = 0; i < size; i++) {
        sum += values.at(i);
    }
    for (unsigned long i = 0; i < size; i++) {
        values[i] = values.at(i) / sum;
    }
}

/**
 * The l2Norm method is used to apply Least Squares, it accumulates second power of each items of values vector
 * and returns the square root of this summation.
 *
 * @return square root of this summation.
 */
double Vector::l2Norm() const{
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += pow(values.at(i), 2);
    }
    return sqrt(sum);
}

/**
 * The cosineSimilarity method takes a Vector v as an input and returns the result of dotProduct(v) / l2Norm() / v.l2Norm().
 * If sizes do not match it throws a VectorSizeMismatch exception.
 *
 * @param v Vector input.
 * @return dotProduct(v) / l2Norm() / v.l2Norm().
 */
double Vector::cosineSimilarity(const Vector& v) const{
    if (size != v.size) {
        throw VectorSizeMismatch();
    }
    return dotProduct(v) / l2Norm() / v.l2Norm();
}

/**
 * Getter for the item at given index of values vector.
 *
 * @param index used to get an item.
 * @return the item at given index.
 */
double Vector::getValue(unsigned long index) const{
    return values.at(index);
}

/**
 * Setter for the setting the value at given index of values vector.
 *
 * @param index to set.
 * @param value is used to set the given index
 */
void Vector::setValue(unsigned long index, double value) {
    values[index] = value;
}

/**
 * The addValue method adds the given value to the item at given index of values vector.
 *
 * @param index to add the given value.
 * @param value value to add to given index.
 */
void Vector::addValue(unsigned long index, double value) {
    values[index] += value;
}

/**
 * The size method returns the size of the values vector.
 *
 * @return size of the values vector.
 */
unsigned long Vector::getSize() const {
    return size;
}

/**
 * The sum method returns the sum of the values vector.
 *
 * @return sum of the values vector.
 */
double Vector::sum() const{
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += values.at(i);
    }
    return sum;
}

/**
 * Swaps two elements with the given indexes.
 * @param index1 First index
 * @param index2 Second index
 */
void Vector::swap(int index1, int index2) {
    double tmp;
    tmp = values.at(index1);
    values[index1] = values.at(index2);
    values[index2] = tmp;
}

/**
 * Writes the vector to an output file.
 * @param outputFile Output file.
 */
void Vector::serialize(ostream &outputFile) {
    outputFile << values.size() << " ";
    for (double value : values){
        outputFile << value << " ";
    }
    outputFile << "\n";
}

/**
 * Reads the vector from an input file.
 * @param inputFile Input file.
 */
Vector::Vector(ifstream &inputFile) {
    inputFile >> size;
    for (int i = 0; i < size; i++){
        double value;
        inputFile >> value;
        values.push_back(value);
    }
}
