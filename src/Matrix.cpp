//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#include <fstream>
#include <random>
#include <complex>
#include "Matrix.h"
#include "MatrixDimensionMismatch.h"
#include "MatrixColumnMismatch.h"
#include "MatrixRowMismatch.h"
#include "MatrixRowColumnMismatch.h"
#include "DeterminantZero.h"
#include "MatrixNotPositiveDefinite.h"
#include "MatrixNotSymmetric.h"
#include "MatrixNotSquare.h"

using namespace std;

/**
 * A constructor of {@link Matrix} class which takes a filename as an input and reads numbers into values {@link array}
 * and row and column variables.
 *
 * @param filename is used to read file.
 */
Matrix::Matrix(const string& fileName) {
    double value;
    ifstream inputStream(fileName, ios::in);
    inputStream >> row;
    inputStream >> col;
    values = new double*[row];
    for (int i = 0; i < row; i++) {
        values[i] = new double[col];
        for (int j = 0; j < col; j++) {
            inputStream >> value;
            values[i][j] = value;
        }
    }
    inputStream.close();
}

Matrix::Matrix(ifstream &inputStream) {
    double value;
    inputStream >> row;
    inputStream >> col;
    values = new double*[row];
    for (int i = 0; i < row; i++) {
        values[i] = new double[col];
        for (int j = 0; j < col; j++) {
            inputStream >> value;
            values[i][j] = value;
        }
    }
}

/**
 * Another constructor of {@link Matrix} class which takes row and column numbers as inputs and creates new values
 * {@link array} with given parameters.
 *
 * @param row is used to create matrix.
 * @param col is used to create matrix.
 */
Matrix::Matrix(int row, int col) {
    values = new double*[row];
    for (int i = 0; i < row; i++){
        values[i] = new double[col];
        for (int j = 0; j < col; j++){
            values[i][j] = 0;
        }
    }
    this->row = row;
    this->col = col;
}

/**
 * Another constructor of {@link Matrix} class which takes row, column, minimum and maximum values as inputs.
 * First it creates new values {@link array} with given row and column numbers. Then fills in the
 * positions with random numbers using minimum and maximum inputs.
 *
 * @param row is used to create matrix.
 * @param col is used to create matrix.
 * @param min minimum value.
 * @param max maximum value.
 */
Matrix::Matrix(int row, int col, double min, double max, default_random_engine randomEngine) {
    values = new double*[row];
    this->row = row;
    this->col = col;
    uniform_real_distribution <> distribution (min, max);
    for (int i = 0; i < row; i++){
        values[i] = new double[col];
        for (int j = 0; j < col; j++){
            values[i][j] = distribution(randomEngine);
        }
    }
}

/**
 * Another constructor of {@link Matrix} class which takes size as input and creates new values {@link array}
 * with using size input and assigns 1 to each element at the diagonal.
 *
 * @param size is used declaring the size of the array.
 */
Matrix::Matrix(int size) {
    row = size;
    col = size;
    values = new double*[row];
    for (int i = 0; i < row; i++){
        values[i] = new double[col];
        for (int j = 0; j < col; j++){
            values[i][j] = 0;
        }
    }
    for (int i = 0; i < row; i++) {
        values[i][i] = 1;
    }
}

/**
 * Another constructor takes two {@link Vector}s v1 and v2 as an input and creates new {@link Matrix} m of [size x size of input v].
 * It loops through the the both values {@link Vector} and given vector's values {@link Vector}, then multiply
 * each item with other with other items and puts to the new {@link Matrix} m.
 *
 * @param v1 Vector input.
 * @param v2 Vector input.
 * @return Matrix that is the multiplication of two vectors.
 */
Matrix::Matrix(const Vector& v1, const Vector& v2) {
    row = v1.getSize();
    col = v2.getSize();
    values = new double*[row];
    for (int i = 0; i < row; i++){
        values[i] = new double[col];
    }
    for (int i = 0; i < v1.getSize(); i++) {
        for (int j = 0; j < v2.getSize(); j++) {
            values[i][j] = v1.getValue(i) * v2.getValue(j);
        }
    }
}

/**
 * The printToFile method takes a fileName as an input and prints values {@link array} into the file.
 *
 * @param fileName String input to write to file.
 */
void Matrix::printToFile(const string& fileName) const{
    ofstream outputStream(fileName, ios::out);
    for (int i = 0; i < row; i++) {
        outputStream << values[i][0];
        for (int j = 1; j < col; j++) {
            outputStream << values[i][j];
        }
        outputStream << "\n";
    }
    outputStream.close();
}

/**
 * The getter for the index at given rowNo and colNo of values {@link array}.
 *
 * @param rowNo integer input for row number.
 * @param colNo integer input for column number.
 * @return item at given index of values {@link array}.
 */
double Matrix::getValue(int rowNo, int colNo) const{
    return values[rowNo][colNo];
}

/**
 * The setter for the value at given index of values {@link aArray}.
 *
 * @param rowNo integer input for row number.
 * @param colNo integer input for column number.
 * @param value is used to set at given index.
 */
void Matrix::setValue(int rowNo, int colNo, double value) {
    values[rowNo][colNo] = value;
}

/**
 * The addValue method adds the given value to the item at given index of values {@link java.lang.reflect.Array}.
 *
 * @param rowNo integer input for row number.
 * @param colNo integer input for column number.
 * @param value is used to add to given item at given index.
 */
void Matrix::addValue(int rowNo, int colNo, double value) {
    values[rowNo][colNo] += value;
}

/**
 * The increment method adds 1 to the item at given index of values {@link array}.
 *
 * @param rowNo integer input for row number.
 * @param colNo integer input for column number.
 */
void Matrix::increment(int rowNo, int colNo) {
    values[rowNo][colNo] += 1;
}

/**
 * The getter for the row variable.
 *
 * @return row number.
 */
int Matrix::getRow() const{
    return row;
}

/**
 * The getRow method returns the vector of values {@link array} at given row input.
 *
 * @param _row integer input for row number.
 * @return Vector of values {@link array} at given row input.
 */
Vector Matrix::getRow(int _row) const{
    return Vector(values[_row], col);
}

/**
 * The getColumn method creates an {@link vector} and adds items at given column number of values {@link array}
 * to the {@link vector}.
 *
 * @param column integer input for column number.
 * @return Vector of given column number.
 */
vector<double> Matrix::getColumn(int column) const{
    vector<double> vector;
    vector.reserve(row);
    for (int i = 0; i < row; i++) {
        vector.push_back(values[i][column]);
    }
    return vector;
}

/**
 * The getter for column variable.
 *
 * @return column variable.
 */
int Matrix::getColumn() const{
    return col;
}

/**
 * The columnWiseNormalize method, first accumulates items column by column then divides items by the summation.
 */
void Matrix::columnWiseNormalize() {
    for (int i = 0; i < row; i++) {
        double sum = 0.0;
        for (int j = 0; j < col; j++) {
            sum += values[i][j];
        }
        for (int j = 0; j < col; j++) {
            values[i][j] /= sum;
        }
    }
}

/**
 * The multiplyWithConstant method takes a constant as an input and multiplies each item of values {@link array}
 * with given constant.
 *
 * @param constant value to multiply items of values {@link array}.
 */
void Matrix::multiplyWithConstant(double constant) {
    int i;
    for (i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            values[i][j] *= constant;
        }
    }
}

/**
 * The divideByConstant method takes a constant as an input and divides each item of values {@link array}
 * with given constant.
 *
 * @param constant value to divide items of values {@link array}.
 */
void Matrix::divideByConstant(double constant) {
    int i;
    for (i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            values[i][j] /= constant;
        }
    }
}

/**
 * The add method takes a {@link Matrix} as an input and accumulates values {@link array} with the
 * corresponding items of given Matrix. If the sizes of both Matrix and values {@link array} do not match,
 * it throws {@link MatrixDimensionMismatch} exception.
 *
 * @param m Matrix type input.
 */
void Matrix::add(const Matrix& m) {
    int i;
    if (row != m.row || col != m.col) {
        throw MatrixDimensionMismatch();
    }
    for (i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            values[i][j] += m.values[i][j];
        }
    }
}

/**
 * The add method which takes a row number and a Vector as inputs. It sums up the corresponding values at the given row of
 * values {@link array} and given {@link Vector}. If the sizes of both Matrix and values
 * {@link array} do not match, it throws {@link MatrixColumnMismatch} exception.
 *
 * @param rowNo integer input for row number.
 * @param v     Vector type input.
 */
void Matrix::add(int rowNo, const Vector& v) {
    if (col != v.getSize()) {
        throw MatrixColumnMismatch();
    }
    for (int i = 0; i < col; i++){
        values[rowNo][i] += v.getValue(i);
    }
}

/**
 * The subtract method takes a {@link Matrix} as an input and subtracts from values {@link array} the
 * corresponding items of given Matrix. If the sizes of both Matrix and values {@link aArray} do not match,
 * it throws {@link MatrixDimensionMismatch} exception.
 *
 * @param m Matrix type input.
 */
void Matrix::subtract(const Matrix& m) {
    int i;
    if (row != m.row || col != m.col) {
        throw MatrixDimensionMismatch();
    }
    for (i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            values[i][j] -= m.values[i][j];
        }
    }
}

/**
 * The multiplyWithVectorFromLeft method takes a Vector as an input and creates a result {@link array}.
 * Then, multiplies values of input Vector starting from the left side with the values {@link array},
 * accumulates the multiplication, and assigns to the result {@link array}. If the sizes of both Vector
 * and row number do not match, it throws {@link MatrixRowMismatch} exception.
 *
 * @param v {@link Vector} type input.
 * @return Vector that holds the result.
 */
Vector Matrix::multiplyWithVectorFromLeft(const Vector& v) const{
    if (row != v.getSize()) {
        throw MatrixRowMismatch();
    }
    vector<double> result((unsigned long) col);
    for (int i = 0; i < col; i++) {
        result[i] = 0.0;
        for (unsigned long j = 0; j < row; j++) {
            result[i] += v.getValue(j) * values[j][i];
        }
    }
    return Vector(result);
}

/**
 * The multiplyWithVectorFromRight method takes a Vector as an input and creates a result {@link array}.
 * Then, multiplies values of input Vector starting from the right side with the values {@link array},
 * accumulates the multiplication, and assigns to the result {@link array}. If the sizes of both Vector
 * and row number do not match, it throws {@link MatrixColumnMismatch} exception.
 *
 * @param v {@link Vector} type input.
 * @return Vector that holds the result.
 */
Vector Matrix::multiplyWithVectorFromRight(const Vector& v) const{
    if (col != v.getSize()) {
        throw new MatrixColumnMismatch();
    }
    vector<double> result((unsigned long) row);
    for (int i = 0; i < row; i++) {
        result[i] = 0;
        for (int j = 0; j < col; j++){
            result[i] += values[i][j] * v.getValue(j);
        }
    }
    return Vector(result);
}

/**
 * The columnSum method takes a column number as an input and accumulates items at given column number of values
 * {@link array}.
 *
 * @param columnNo Column number input.
 * @return summation of given column of values {@link array}.
 */
double Matrix::columnSum(int columnNo) const{
    double sum = 0;
    for (int i = 0; i < row; i++) {
        sum += values[i][columnNo];
    }
    return sum;
}

/**
 * The sumOfRows method creates a mew result {@link Vector} and adds the result of columnDum method's corresponding
 * index to the newly created result {@link Vector}.
 *
 * @return Vector that holds column sum.
 */
Vector Matrix::sumOfRows() const{
    Vector result = Vector(0, 0.0);
    for (int i = 0; i < col; i++) {
        result.add(columnSum(i));
    }
    return result;
}

/**
 * The rowSum method takes a row number as an input and accumulates items at given row number of values
 * {@link array}.
 *
 * @param rowNo Row number input.
 * @return summation of given row of values {@link array}.
 */
double Matrix::rowSum(int rowNo) const{
    double sum = 0;
    for (int i = 0; i < col; i++){
        sum += values[rowNo][i];
    }
    return sum;
}

/**
 * The multiply method takes a {@link Matrix} as an input. First it creates a result {@link Matrix} and puts the
 * accumulatated multiplication of values {@link array} and given {@link Matrix} into result
 * {@link Matrix}. If the size of Matrix's row size and values {@link array}'s column size do not match,
 * it throws {@link MatrixRowColumnMismatch} exception.
 *
 * @param m Matrix type input.
 * @return result {@link Matrix}.
 */
Matrix Matrix::multiply(const Matrix& m) const{
    int i, j, k;
    double sum;
    Matrix result(row, m.col);
    if (col != m.row) {
        throw MatrixRowColumnMismatch();
    }
    for (i = 0; i < row; i++) {
        for (j = 0; j < m.col; j++) {
            sum = 0.0;
            for (k = 0; k < col; k++) {
                sum += values[i][k] * m.values[k][j];
            }
            result.values[i][j] = sum;
        }
    }
    return result;
}

/**
 * The elementProduct method takes a {@link Matrix} as an input and performs element wise multiplication. Puts result
 * to the newly created Matrix. If the size of Matrix's row and column size does not match with the values
 * {@link array}'s row and column size, it throws {@link MatrixDimensionMismatch} exception.
 *
 * @param m Matrix type input.
 * @return result {@link Matrix}.
 */
Matrix Matrix::elementProduct(const Matrix& m) const{
    int i;
    if (row != m.row || col != m.col) {
        throw MatrixDimensionMismatch();
    }
    Matrix result(row, m.col);
    for (i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            result.values[i][j] = values[i][j] * m.values[i][j];
        }
    }
    return result;
}

/**
 * The sumOfElements method accumulates all the items in values {@link array} and
 * returns this summation.
 *
 * @return sum of the items of values {@link array}.
 */
double Matrix::sumOfElements() const{
    int i;
    double sum = 0.0;
    for (i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            sum += values[i][j];
        }
    }
    return sum;
}

/**
 * The trace method accumulates items of values {@link java.lang.reflect.Array} at the diagonal.
 *
 * @return sum of items at diagonal.
 */
double Matrix::trace() const{
    int i;
    if (row != col){
        throw MatrixNotSquare();
    }
    double sum = 0.0;
    for (i = 0; i < row; i++) {
        sum += values[i][i];
    }
    return sum;
}

/**
 * The transpose method creates a new {@link Matrix}, then takes the transpose of values {@link array}
 * and puts transposition to the {@link Matrix}.
 *
 * @return Matrix type output.
 */
Matrix Matrix::transpose() {
    int i, j;
    Matrix result(col, row);
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            result.values[j][i] = values[i][j];
        }
    }
    return result;
}

/**
 * The partial method takes 4 integer inputs; rowstart, rowend, colstart, colend and creates a {@link Matrix} size of
 * rowend - rowstart + 1 x colend - colstart + 1. Then, puts corresponding items of values {@link array}
 * to the new result {@link Matrix}.
 *
 * @param rowstart integer input for defining starting index of row.
 * @param rowend   integer input for defining ending index of row.
 * @param colstart integer input for defining starting index of column.
 * @param colend   integer input for defining ending index of column.
 * @return result Matrix.
 */
Matrix Matrix::partial(int rowstart, int rowend, int colstart, int colend) {
    int i, j;
    Matrix result(rowend - rowstart + 1, colend - colstart + 1);
    for (i = rowstart; i <= rowend; i++)
        for (j = colstart; j <= colend; j++)
            result.values[i - rowstart][j - colstart] = values[i][j];
    return result;
}

/**
 * The isSymmetric method compares each item of values {@link array} at positions (i, j) with (j, i)
 * and returns true if they are equal, false otherwise.
 *
 * @return true if items are equal, false otherwise.
 */
bool Matrix::isSymmetric() {
    if (row != col){
        throw MatrixNotSquare();
    }
    for (int i = 0; i < row - 1; i++) {
        for (int j = i + 1; j < row; j++) {
            if (values[i][j] != values[j][i]) {
                return false;
            }
        }
    }
    return true;
}

/**
 * The determinant method first creates a new {@link array}, and copies the items of  values
 * {@link array} into new {@link array}. Then, calculates the determinant of this
 * new {@link array}.
 *
 * @return determinant of values {@link array}.
 */
double Matrix::determinant() {
    if (row != col){
        throw MatrixNotSquare();
    }
    int i, j, k;
    double ratio, det = 1.0;
    Matrix copy = clone();
    for (i = 0; i < row; i++) {
        det *= copy.values[i][i];
        if (det == 0.0)
            break;
        for (j = i + 1; j < row; j++) {
            ratio = copy.values[j][i] / copy.values[i][i];
            for (k = i; k < col; k++)
                copy.values[j][k] -= copy.values[i][k] * ratio;
        }
    }
    return det;
}

/**
 * The inverse method finds the inverse of values {@link array}.
 */
void Matrix::inverse() {
    if (row != col){
        throw MatrixNotSquare();
    }
    double big;
    double dum, pivinv;
    int i, icol, irow, k, l, ll;
    Matrix b(row);
    int* indxc, *indxr, *ipiv;
    indxc = new int[row];
    indxr = new int[row];
    ipiv = new int[row];
    for (int j = 0; j < row; j++)
        ipiv[j] = 0;
    for (i = 1; i <= row; i++) {
        big = 0.0;
        irow = -1;
        icol = -1;
        for (int j = 1; j <= row; j++)
            if (ipiv[j - 1] != 1)
                for (k = 1; k <= row; k++)
                    if (ipiv[k - 1] == 0)
                        if (abs(values[j - 1][k - 1]) >= big) {
                            big = abs(values[j - 1][k - 1]);
                            irow = j;
                            icol = k;
                        }
        if (irow == -1 || icol == -1)
            throw DeterminantZero();
        ipiv[icol - 1] = ipiv[icol - 1] + 1;
        if (irow != icol) {
            auto* dummy = new double[col];
            memcpy(dummy, values[irow - 1], col * sizeof(double));
            memcpy(values[irow - 1], values[icol - 1], col * sizeof(double));
            memcpy(values[icol - 1], dummy, col * sizeof(double));
            memcpy(dummy, b.values[irow - 1], col * sizeof(double));
            memcpy(b.values[irow - 1], b.values[icol - 1], col * sizeof(double));
            memcpy(b.values[icol - 1], dummy, col * sizeof(double));
            delete[] dummy;
        }
        indxr[i - 1] = irow;
        indxc[i - 1] = icol;
        if (values[icol - 1][icol - 1] == 0)
            throw DeterminantZero();
        pivinv = (1.0) / (values[icol - 1][icol - 1]);
        values[icol - 1][icol - 1] = 1.0;
        for (int j = 0; j < col; j++){
            values[icol - 1][j] *= pivinv;
            b.values[icol - 1][j] *= pivinv;
        }
        for (ll = 1; ll <= row; ll++)
            if (ll != icol) {
                dum = values[ll - 1][icol - 1];
                values[ll - 1][icol - 1] = 0.0;
                for (l = 1; l <= row; l++)
                    values[ll - 1][l - 1] -= values[icol - 1][l - 1] * dum;
                for (l = 1; l <= row; l++)
                    b.values[ll - 1][l - 1] -= b.values[icol - 1][l - 1] * dum;
            }
    }
    for (l = row; l >= 1; l--){
        if (indxr[l - 1] != indxc[l - 1]){
            for (k = 1; k <= row; k++) {
                double tmp = values[k - 1][indxr[l - 1] - 1];
                values[k - 1][indxr[l - 1] - 1] = values[k - 1][indxc[l - 1] - 1];
                values[k - 1][indxc[l - 1] - 1] = tmp;
            }
        }
    }
    delete[] indxc;
    delete[] indxr;
    delete[] ipiv;
}

/**
 * The choleskyDecomposition method creates a new {@link Matrix} and puts the Cholesky Decomposition of values Array
 * into this {@link Matrix}. Also, it throws {@link MatrixNotSymmetric} exception if it is not symmetric and
 * {@link MatrixNotPositiveDefinite} exception if the summation is negative.
 *
 * @return Matrix type output.
 */
Matrix Matrix::choleskyDecomposition() {
    int i, j, k;
    double sum;
    if (!isSymmetric()) {
        throw MatrixNotSymmetric();
    }
    Matrix b(row, col);
    for (i = 0; i < row; i++) {
        for (j = i; j < row; j++) {
            sum = values[i][j];
            for (k = i - 1; k >= 0; k--)
                sum -= values[i][k] * values[j][k];
            if (i == j) {
                if (sum <= 0.0)
                    throw MatrixNotPositiveDefinite();
                b.values[i][i] = sqrt(sum);
            } else
                b.values[j][i] = sum / b.values[i][i];
        }
    }
    return b;
}

/**
 * The rotate method rotates values {@link array} according to given inputs.
 *
 * @param s   double input.
 * @param tau double input.
 * @param i   integer input.
 * @param j   integer input.
 * @param k   integer input.
 * @param l   integer input.
 */
void Matrix::rotate(double s, double tau, int i, int j, int k, int l) {
    double g = values[i][j];
    double h = values[k][l];
    values[i][j] = g - s * (h + g * tau);
    values[k][l] = h + s * (g - h * tau);
}

bool comparator(const Eigenvector& i, const Eigenvector& j) { return (i.getEigenValue() < j.getEigenValue()); }

/**
 * The characteristics method finds and returns a sorted {@link vector} of {@link Eigenvector}s. And it throws
 * {@link MatrixNotSymmetric} exception if it is not symmetric.
 *
 * @return a sorted {@link vector} of {@link Eigenvector}s.
 */
vector<Eigenvector> Matrix::characteristics() {
    int j, iq, ip, i;
    double threshold, theta, tau, t, sm, s, h, g, c;
    if (!isSymmetric()) {
        throw MatrixNotSymmetric();
    }
    Matrix matrix1 = this->clone();
    Matrix v(row);
    auto* d = new double[row];
    auto* b = new double[row];
    auto* z = new double[row];
    double EPS = 0.000000000000000001;
    for (ip = 0; ip < row; ip++) {
        b[ip] = d[ip] = matrix1.values[ip][ip];
        z[ip] = 0.0;
    }
    for (i = 1; i <= 50; i++) {
        sm = 0.0;
        for (ip = 0; ip < row - 1; ip++)
            for (iq = ip + 1; iq < row; iq++)
                sm += abs(matrix1.values[ip][iq]);
        if (sm == 0.0) {
            break;
        }
        if (i < 4)
            threshold = 0.2 * sm / pow(row, 2);
        else
            threshold = 0.0;
        for (ip = 0; ip < row - 1; ip++) {
            for (iq = ip + 1; iq < row; iq++) {
                g = 100.0 * abs(matrix1.values[ip][iq]);
                if (i > 4 && g <= EPS * abs(d[ip]) && g <= EPS * abs(d[iq])) {
                    matrix1.values[ip][iq] = 0.0;
                } else {
                    if (abs(matrix1.values[ip][iq]) > threshold) {
                        h = d[iq] - d[ip];
                        if (g <= EPS * abs(h)) {
                            t = matrix1.values[ip][iq] / h;
                        } else {
                            theta = 0.5 * h / matrix1.values[ip][iq];
                            t = 1.0 / (abs(theta) + sqrt(1.0 + pow(theta, 2)));
                            if (theta < 0.0) {
                                t = -t;
                            }
                        }
                        c = 1.0 / sqrt(1 + pow(t, 2));
                        s = t * c;
                        tau = s / (1.0 + c);
                        h = t * matrix1.values[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        matrix1.values[ip][iq] = 0.0;
                        for (j = 0; j < ip; j++) {
                            matrix1.rotate(s, tau, j, ip, j, iq);
                        }
                        for (j = ip + 1; j < iq; j++) {
                            matrix1.rotate(s, tau, ip, j, j, iq);
                        }
                        for (j = iq + 1; j < row; j++) {
                            matrix1.rotate(s, tau, ip, j, iq, j);
                        }
                        for (j = 0; j < row; j++) {
                            v.rotate(s, tau, j, ip, j, iq);
                        }
                    }
                }
            }
        }
        for (ip = 0; ip < row; ip++) {
            b[ip] = b[ip] + z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    vector<Eigenvector> result;
    for (i = 0; i < row; i++) {
        if (d[i] > 0) {
            result.emplace_back(d[i], v.getColumn(i));
        }
    }
    sort(result.begin(), result.end(), comparator);
    delete[] d;
    delete[] b;
    delete[] z;
    return result;
}

Matrix Matrix::clone() {
    Matrix result = Matrix(row, col);
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
            result.values[i][j] = values[i][j];
        }
    }
    return result;
}

Matrix Matrix::sum(const Matrix &m) const{
    int i, j;
    Matrix result = Matrix(row, col);
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            result.values[i][j] = values[i][j] + m.values[i][j];
        }
    }
    return result;
}

Matrix Matrix::difference(const Matrix &m) const {
    int i, j;
    Matrix result = Matrix(row, col);
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            result.values[i][j] = values[i][j] - m.values[i][j];
        }
    }
    return result;
}
