//
// Created by Olcay Taner Yıldız on 26.11.2018.
//

#include "VectorSizeMismatch.h"

/**
 * The overridden toString method returns 'Number of items in both vectors must be the same' String.
 *
 * @return 'Number of items in both vectors must be the same' String.
 */
const char *VectorSizeMismatch::what() const noexcept {
    return "Number of items in both vectors must be the same";
}
