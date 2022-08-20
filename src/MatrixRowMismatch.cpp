//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "MatrixRowMismatch.h"

/**
 * The overridden toString method returns 'Number of rows of the matrix should be equal to the size of the vector' String.
 *
 * @return 'Number of rows of the matrix should be equal to the size of the vector' String.
 */
const char *MatrixRowMismatch::what() const noexcept {
    return "Number of rows of the matrix should be equal to the size of the vector.";
}
