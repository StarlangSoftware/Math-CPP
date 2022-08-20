//
// Created by Olcay Taner Yıldız on 4.12.2018.
//

#include "MatrixNotSquare.h"

/**
 * The overridden toString method returns 'Matrix should be square' String.
 *
 * @return 'Matrix should be square' String.
 */
const char *MatrixNotSquare::what() const noexcept {
    return "Matrix should be square.";
}
