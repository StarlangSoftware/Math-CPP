//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "MatrixNotPositiveDefinite.h"

/**
 * The overridden toString method returns 'Matrix should be positive definite' String.
 *
 * @return 'Matrix should be positive definite' String.
 */
const char *MatrixNotPositiveDefinite::what() const noexcept {
    return "Matrix should be positive definite.";
}
