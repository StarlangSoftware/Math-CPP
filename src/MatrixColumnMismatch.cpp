//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "MatrixColumnMismatch.h"

/**
 * The overridden toString method returns 'Number of columns of the matrix should be equal to the size of the vector' String.
 *
 * @return 'Number of columns of the matrix should be equal to the size of the vector' String.
 */
const char *MatrixColumnMismatch::what() const noexcept {
    return "Number of columns of the matrix should be equal to the size of the vector.";
}
