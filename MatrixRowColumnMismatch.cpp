//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "MatrixRowColumnMismatch.h"

/**
 * The overridden toString method returns 'The number of columns of the first matrix should be equal to the number
 * of rows of the second matrix' String.
 *
 * @return 'The number of columns of the first matrix should be equal to the number of rows of the second matrix' String.
 */
const char *MatrixRowColumnMismatch::what() const noexcept {
    return "The number of columns of the first matrix should be equal to the number of rows of the second matrix.";
}
