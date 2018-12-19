//
// Created by Olcay Taner Yıldız on 19.12.2018.
//

#include "Matrix.h"

int main(){
    Matrix matrix(2, 2);
    matrix.increment(0, 0);
    matrix.increment(0, 1);
    matrix.increment(1, 0);
    matrix.increment(1, 1);
    matrix.increment(0, 0);
    matrix.increment(0, 1);
    matrix.increment(1, 0);
    matrix.increment(1, 1);
}