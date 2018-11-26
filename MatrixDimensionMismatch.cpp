//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "MatrixDimensionMismatch.h"

/**
 * The overridden toString method returns 'The number of rows and columns of the first matrix should be equal to the
 * number of rows and columns of the second matrix respectively' String.
 *
 * @return 'The number of rows and columns of the first matrix should be equal to the number of rows and columns of
 * the second matrix respectively' String.
 */
const char *MatrixDimensionMismatch::what() const noexcept {
    return "The number of rows and columns of the first matrix should be equal to the number of rows and columns of the second matrix respectively.";
}
