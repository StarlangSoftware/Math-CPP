//
// Created by Yiğit Demirşan on 23.02.2025.
//

#include "catch.hpp"
#include "../src/Tensor.h"

TEST_CASE("Tensor Initialization from Nested Data") {
    vector<vector<vector<float> > > nested_data = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
    Tensor tensor(nested_data);

    REQUIRE(tensor.get({0, 0, 0}) == 1);
    REQUIRE(tensor.get({1, 1, 1}) == 8);
}

TEST_CASE("Tensor Initialization from Flat Data") {
    vector<float> flat_data = {1, 2, 3, 4, 5, 6};
    vector<int> shape = {2, 3};
    Tensor tensor(flat_data, shape);

    REQUIRE(tensor.get({0, 0}) == 1);
    REQUIRE(tensor.get({1, 2}) == 6);
}

TEST_CASE("Get and Set Methods") {
    vector<float> flat_data = {1, 2, 3, 4, 5, 6};
    vector<int> shape = {2, 3};
    Tensor tensor(flat_data, shape);

    tensor.set({0, 1}, 10);
    REQUIRE(tensor.get({0, 1}) == 10);
}

TEST_CASE("Reshape Tensor") {
    vector<float> flat_data = {1, 2, 3, 4};
    vector<int> shape = {2, 2};
    Tensor tensor(flat_data, shape);

    Tensor reshaped = tensor.reshape({4});
    REQUIRE(reshaped.get({2}) == 3);
}

TEST_CASE("Transpose Tensor") {
    vector<float> flat_data = {1, 2, 3, 4, 5, 6};
    vector<int> shape = {2, 3};
    Tensor tensor(flat_data, shape);

    Tensor transposed = tensor.transpose({1, 0});
    REQUIRE(transposed.get({0, 1}) == 3);
    REQUIRE(transposed.get({2, 1}) == 6);
}

TEST_CASE("Addition Operator") {
    vector<float> flat_data1 = {1, 2, 3, 4};
    vector<float> flat_data2 = {5, 6, 7, 8};
    vector<int> shape = {2, 2};

    Tensor tensor1(flat_data1, shape);
    Tensor tensor2(flat_data2, shape);

    Tensor result = tensor1 + tensor2;

    REQUIRE(result.get({0, 0}) == 6);
    REQUIRE(result.get({1, 1}) == 12);
}

TEST_CASE("Subtraction Operator") {
    vector<float> flat_data1 = {5, 6, 7, 8};
    vector<float> flat_data2 = {1, 2, 3, 4};
    vector<int> shape = {2, 2};

    Tensor tensor1(flat_data1, shape);
    Tensor tensor2(flat_data2, shape);

    Tensor result = tensor1 - tensor2;

    REQUIRE(result.get({0, 0}) == 4);
    REQUIRE(result.get({1, 1}) == 4);
}

TEST_CASE("Multiplication Operator") {
    vector<float> flat_data1 = {1, 2, 3, 4};
    vector<float> flat_data2 = {2, 3, 4, 5};
    vector<int> shape = {2, 2};

    Tensor tensor1(flat_data1, shape);
    Tensor tensor2(flat_data2, shape);

    Tensor result = tensor1 * tensor2;

    REQUIRE(result.get({0, 0}) == 2);
    REQUIRE(result.get({1, 1}) == 20);
}

TEST_CASE("Dot Product") {
    vector<float> flat_data1 = {1, 2, 3, 4};
    vector<float> flat_data2 = {1, 2, 3, 4};
    vector<int> shape1 = {2, 2};
    vector<int> shape2 = {2, 2};

    Tensor tensor1(flat_data1, shape1);
    Tensor tensor2(flat_data2, shape2);

    Tensor result = tensor1.dot(tensor2);

    REQUIRE(result.get({0, 0}) == 7);
    REQUIRE(result.get({1, 1}) == 25);
}

TEST_CASE("ToString Method") {
    vector<float> flat_data = {1, 2, 3, 4};
    vector<int> shape = {2, 2};
    Tensor tensor(flat_data, shape);

    string tensor_str = tensor.to_string();
    REQUIRE(tensor_str.find("Tensor(shape=[2, 2], data=[1, 2, 3, 4])") != string::npos);
}

TEST_CASE("Invalid Reshape") {
    vector<float> flat_data = {1, 2, 3, 4};
    vector<int> shape = {2, 2};
    Tensor tensor(flat_data, shape);

    REQUIRE_THROWS_AS(tensor.reshape({3, 2}), invalid_argument);
}

TEST_CASE("Invalid Dot Product") {
    vector<float> flat_data1 = {1, 2, 3, 4};
    vector<float> flat_data2 = {1, 2, 3, 4};
    vector<int> shape1 = {2, 2};
    vector<int> shape2 = {3, 2};

    Tensor tensor1(flat_data1, shape1);
    Tensor tensor2(flat_data2, shape2);

    REQUIRE_THROWS_AS(tensor1.dot(tensor2), invalid_argument);
}
