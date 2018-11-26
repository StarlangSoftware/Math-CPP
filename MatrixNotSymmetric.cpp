//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "MatrixNotSymmetric.h"

/**
 * The overridden toString method returns 'Matrix should be symmetric' String.
 *
 * @return 'Matrix should be symmetric' String.
 */
const char *MatrixNotSymmetric::what() const noexcept {
    return "Matrix should be symmetric.";
}
