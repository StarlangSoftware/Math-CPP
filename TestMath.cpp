//
// Created by Olcay Taner Yıldız on 19.12.2018.
//

#include <fstream>
#include "Matrix.h"

int main(){
    Vector v1(100, 2);
    ofstream output;
    ifstream input;
    output.open("deneme.txt", ofstream::out);
    v1.serialize(output);
    output.close();
    input.open("deneme.txt", ifstream::in);
    Vector v2(input);
    input.close();
    Matrix m1(10, 10, 0.0, 1.0);
    m1.inverse();
    output.open("deneme.txt", ofstream::out);
    m1.serialize(output);
    output.close();
    input.open("deneme.txt", ifstream::in);
    Matrix m2(input);
    input.close();
}