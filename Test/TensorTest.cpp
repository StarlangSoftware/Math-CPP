//
// Created by Yiğit Demirşan on 23.02.2025.
//

#include "catch.hpp"
#include "../src/Tensor.h"

using namespace std;

/**
 * Checks whether the shape of the tensor is equal to the expected shape.
 */
static bool shapeEquals(const Tensor& tensor, const initializer_list<int> expected) {
    if (tensor.getShapeSize() != static_cast<int>(expected.size())) {
        return false;
    }
    int i = 0;
    for (const int dimension : expected) {
        if (tensor.getShape()[i++] != dimension) {
            return false;
        }
    }
    return true;
}

/**
 * Checks whether the data of the tensor is equal to the expected data.
 */
static bool dataEquals(const Tensor& tensor, const initializer_list<double> expected) {
    if (tensor.getSize() != static_cast<int>(expected.size())) {
        return false;
    }
    int i = 0;
    for (const double value : expected) {
        if (abs(tensor.getData()[i++] - value) > 1e-9) {
            return false;
        }
    }
    return true;
}

TEST_CASE("TensorTest-testConstructorWithInferredShape") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE(shapeEquals(tensor, {2, 2}));
    REQUIRE(dataEquals(tensor, {1.0, 2.0, 3.0, 4.0}));
}

TEST_CASE("TensorTest-testConstructorWithExplicitShape") {
    const Tensor tensor = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    REQUIRE(shapeEquals(tensor, {2, 2}));
    REQUIRE(dataEquals(tensor, {1.0, 2.0, 3.0, 4.0}));
}

TEST_CASE("TensorTest-testConstructorWithDataArray") {
    const auto data = new double[4]{1.0, 2.0, 3.0, 4.0};
    const Tensor tensor = Tensor(data, {2, 2});
    REQUIRE(shapeEquals(tensor, {2, 2}));
    REQUIRE(dataEquals(tensor, {1.0, 2.0, 3.0, 4.0}));
    REQUIRE(tensor.getData() == data);
}

TEST_CASE("TensorTest-testConstructorShapeMismatch") {
    REQUIRE_THROWS_AS(Tensor({1.0, 2.0, 3.0}, {2, 2}), invalid_argument);
}

TEST_CASE("TensorTest-testGetShape") {
    const Tensor tensor = Tensor({1.0, 2.0, 3.0}, {3});
    REQUIRE(shapeEquals(tensor, {3}));
    REQUIRE(tensor.getShapeSize() == 1);
    REQUIRE(tensor.getSize() == 3);
}

TEST_CASE("TensorTest-testGetData") {
    const Tensor tensor = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    REQUIRE(dataEquals(tensor, {1.0, 2.0, 3.0, 4.0}));
}

TEST_CASE("TensorTest-testGetValidIndices") {
    const Tensor tensor = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    REQUIRE(tensor.getValue({0, 0}) == 1.0);
    REQUIRE(tensor.getValue({1, 1}) == 4.0);
    const int indices[2] = {1, 0};
    REQUIRE(tensor.getValue(indices) == 3.0);
}

TEST_CASE("TensorTest-testGetOutOfBoundsIndices") {
    const Tensor tensor = Tensor({1.0, 2.0}, {2});
    REQUIRE_THROWS_AS(tensor.getValue({2}), out_of_range);
    REQUIRE_THROWS_AS(tensor.getValue({0, 0}), out_of_range);
}

TEST_CASE("TensorTest-testSetValidIndices") {
    Tensor tensor = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    tensor.setValue({0, 0}, 5.0);
    REQUIRE(tensor.getValue({0, 0}) == 5.0);
    const int indices[2] = {1, 1};
    tensor.setValue(indices, 6.0);
    REQUIRE(tensor.getValue({1, 1}) == 6.0);
}

TEST_CASE("TensorTest-testSetOutOfBoundsIndices") {
    Tensor tensor = Tensor({1.0, 2.0}, {2});
    REQUIRE_THROWS_AS(tensor.setValue({2}, 5.0), out_of_range);
}

TEST_CASE("TensorTest-testCopyDoesNotShareData") {
    Tensor tensor = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    Tensor copied = tensor;
    copied.setValue({0, 0}, 9.0);
    REQUIRE(tensor.getValue({0, 0}) == 1.0);
    REQUIRE(copied.getValue({0, 0}) == 9.0);
    Tensor assigned = Tensor({0.0}, {1});
    assigned = tensor;
    assigned.setValue({1, 1}, 7.0);
    REQUIRE(tensor.getValue({1, 1}) == 4.0);
    REQUIRE(assigned.getValue({1, 1}) == 7.0);
}

TEST_CASE("TensorTest-testReshapeValid") {
    const Tensor tensor = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    const Tensor reshaped = tensor.reshape({4});
    REQUIRE(shapeEquals(reshaped, {4}));
    REQUIRE(dataEquals(reshaped, {1.0, 2.0, 3.0, 4.0}));
}

TEST_CASE("TensorTest-testReshapeInvalid") {
    const Tensor tensor = Tensor({1.0, 2.0, 3.0}, {3});
    REQUIRE_THROWS_AS(tensor.reshape({2, 2}), invalid_argument);
}

TEST_CASE("TensorTest-testTransposeNoAxes") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    const Tensor transposed = tensor.transpose();
    REQUIRE(shapeEquals(transposed, {2, 2}));
    REQUIRE(dataEquals(transposed, {1.0, 3.0, 2.0, 4.0}));
}

TEST_CASE("TensorTest-testTransposeWithAxes") {
    const Tensor tensor = Tensor({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    const Tensor transposed = tensor.transpose({1, 0});
    REQUIRE(shapeEquals(transposed, {3, 2}));
    REQUIRE(dataEquals(transposed, {1.0, 4.0, 2.0, 5.0, 3.0, 6.0}));
}

TEST_CASE("TensorTest-testTransposeInvalidAxesLength") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE_THROWS_AS(tensor.transpose({0}), invalid_argument);
}

TEST_CASE("TensorTest-testTransposeInvalidAxesValue") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE_THROWS_AS(tensor.transpose({0, 2}), invalid_argument);
}

TEST_CASE("TensorTest-testBroadcastToValid") {
    const Tensor tensor = Tensor({1.0, 2.0}, {1, 2});
    const Tensor broadcasted = tensor.broadcastTo({2, 2});
    REQUIRE(shapeEquals(broadcasted, {2, 2}));
    REQUIRE(dataEquals(broadcasted, {1.0, 2.0, 1.0, 2.0}));
}

TEST_CASE("TensorTest-testBroadcastToInvalid") {
    const Tensor tensor = Tensor({1.0, 2.0}, {2});
    REQUIRE_THROWS_AS(tensor.broadcastTo({3}), invalid_argument);
}

TEST_CASE("TensorTest-testAddSameShape") {
    const Tensor tensor1 = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    const Tensor tensor2 = Tensor({5.0, 6.0, 7.0, 8.0}, {2, 2});
    const Tensor sum = tensor1.add(tensor2);
    REQUIRE(shapeEquals(sum, {2, 2}));
    REQUIRE(dataEquals(sum, {6.0, 8.0, 10.0, 12.0}));
}

TEST_CASE("TensorTest-testAddWithBroadcasting") {
    const Tensor tensor1 = Tensor({1.0, 3.0}, {1, 1, 2});
    const Tensor tensor2 = Tensor({2.0, 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 3.0}, {1, 6, 2});
    const Tensor sum = tensor1.add(tensor2);
    REQUIRE(shapeEquals(sum, {1, 6, 2}));
    REQUIRE(dataEquals(sum, {3.0, 6.0, 3.0, 6.0, 3.0, 6.0, 3.0, 6.0, 3.0, 6.0, 3.0, 6.0}));
}

TEST_CASE("TensorTest-testSubtractSameShape") {
    const Tensor tensor1 = Tensor({5.0, 6.0, 7.0, 8.0}, {2, 2});
    const Tensor tensor2 = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    const Tensor difference = tensor1.subtract(tensor2);
    REQUIRE(shapeEquals(difference, {2, 2}));
    REQUIRE(dataEquals(difference, {4.0, 4.0, 4.0, 4.0}));
}

TEST_CASE("TensorTest-testSubtractWithBroadcasting") {
    const Tensor tensor1 = Tensor({5.0, 5.0}, {2, 1});
    const Tensor tensor2 = Tensor({1.0, 2.0}, {1, 2});
    const Tensor difference = tensor1.subtract(tensor2);
    REQUIRE(shapeEquals(difference, {2, 2}));
    REQUIRE(dataEquals(difference, {4.0, 3.0, 4.0, 3.0}));
}

TEST_CASE("TensorTest-testHadamardProductSameShape") {
    const Tensor tensor1 = Tensor({1.0, 2.0, 3.0, 4.0}, {2, 2});
    const Tensor tensor2 = Tensor({5.0, 6.0, 7.0, 8.0}, {2, 2});
    const Tensor product = tensor1.hadamardProduct(tensor2);
    REQUIRE(shapeEquals(product, {2, 2}));
    REQUIRE(dataEquals(product, {5.0, 12.0, 21.0, 32.0}));
}

TEST_CASE("TensorTest-testHadamardProductWithBroadcasting") {
    const Tensor tensor1 = Tensor({1.0, 2.0}, {1, 2});
    const Tensor tensor2 = Tensor({3.0, 4.0}, {2, 1});
    const Tensor product = tensor1.hadamardProduct(tensor2);
    REQUIRE(shapeEquals(product, {2, 2}));
    REQUIRE(dataEquals(product, {3.0, 6.0, 4.0, 8.0}));
}

TEST_CASE("TensorTest-testMatrixMultiply2D") {
    const Tensor tensor1 = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    const Tensor tensor2 = Tensor({{5.0, 6.0}, {7.0, 8.0}});
    const Tensor result = tensor1.multiply(tensor2);
    REQUIRE(shapeEquals(result, {2, 2}));
    REQUIRE(dataEquals(result, {19.0, 22.0, 43.0, 50.0}));
}

TEST_CASE("TensorTest-testMatrixMultiplyNonSquare") {
    const Tensor tensor1 = Tensor({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    const Tensor tensor2 = Tensor({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}});
    const Tensor result = tensor1.multiply(tensor2);
    REQUIRE(shapeEquals(result, {2, 2}));
    REQUIRE(dataEquals(result, {22.0, 28.0, 49.0, 64.0}));
}

TEST_CASE("TensorTest-testMatrixMultiply3D") {
    const Tensor tensor1 = Tensor({{{1.0, 2.0}, {3.0, 4.0}}, {{9.0, 10.0}, {11.0, 12.0}}});
    const Tensor tensor2 = Tensor({{{5.0, 6.0}, {7.0, 8.0}}, {{13.0, 14.0}, {15.0, 16.0}}});
    const Tensor result = tensor1.multiply(tensor2);
    REQUIRE(shapeEquals(result, {2, 2, 2}));
    REQUIRE(dataEquals(result, {19.0, 22.0, 43.0, 50.0, 267.0, 286.0, 323.0, 346.0}));
}

TEST_CASE("TensorTest-testMatrixMultiplyWithBroadcasting") {
    const Tensor tensor1 = Tensor({{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}});
    const Tensor tensor2 = Tensor({1.0, 0.0, 0.0, 1.0}, {1, 2, 2});
    const Tensor result = tensor1.multiply(tensor2);
    REQUIRE(shapeEquals(result, {2, 2, 2}));
    REQUIRE(dataEquals(result, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}));
}

TEST_CASE("TensorTest-testMatrixMultiplyInvalidShapes") {
    const Tensor tensor1 = Tensor({1.0, 2.0, 3.0}, {3});
    const Tensor tensor2 = Tensor({{4.0, 5.0}, {6.0, 7.0}});
    REQUIRE_THROWS_AS(tensor1.multiply(tensor2), invalid_argument);
}

TEST_CASE("TensorTest-testMatrixMultiplyUnsupportedDimensions") {
    const Tensor tensor1 = Tensor({1.0, 2.0}, {2});
    const Tensor tensor2 = Tensor({3.0, 4.0}, {2});
    REQUIRE_THROWS_AS(tensor1.multiply(tensor2), invalid_argument);
}

TEST_CASE("TensorTest-testConcat") {
    const Tensor tensor1 = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    const Tensor tensor2 = Tensor({{5.0, 6.0}, {7.0, 8.0}});
    Tensor concatenated = tensor1.concat(tensor2, 0);
    REQUIRE(shapeEquals(concatenated, {4, 2}));
    REQUIRE(dataEquals(concatenated, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}));
    concatenated = tensor1.concat(tensor2, 1);
    REQUIRE(shapeEquals(concatenated, {2, 4}));
    REQUIRE(dataEquals(concatenated, {1.0, 2.0, 5.0, 6.0, 3.0, 4.0, 7.0, 8.0}));
}

TEST_CASE("TensorTest-testConcatInvalid") {
    const Tensor tensor1 = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    const Tensor tensor2 = Tensor({1.0, 2.0}, {2});
    const Tensor tensor3 = Tensor({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, {3, 2});
    REQUIRE_THROWS_AS(tensor1.concat(tensor1, 2), invalid_argument);
    REQUIRE_THROWS_AS(tensor1.concat(tensor2, 0), invalid_argument);
    REQUIRE_THROWS_AS(tensor1.concat(tensor3, 1), invalid_argument);
}

TEST_CASE("TensorTest-testGet") {
    const Tensor tensor = Tensor({{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}});
    Tensor subTensor = tensor.get({1});
    REQUIRE(shapeEquals(subTensor, {2, 2}));
    REQUIRE(dataEquals(subTensor, {5.0, 6.0, 7.0, 8.0}));
    subTensor = tensor.get({0, 1});
    REQUIRE(shapeEquals(subTensor, {2}));
    REQUIRE(dataEquals(subTensor, {3.0, 4.0}));
}

TEST_CASE("TensorTest-testGetInvalid") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE_THROWS_AS(tensor.get({0, 1}), invalid_argument);
    REQUIRE_THROWS_AS(tensor.get({2}), out_of_range);
}

TEST_CASE("TensorTest-testPartialValid") {
    const Tensor tensor = Tensor({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}});
    const Tensor partial = tensor.partial({0, 1}, {2, 3});
    REQUIRE(shapeEquals(partial, {2, 2}));
    REQUIRE(dataEquals(partial, {2.0, 3.0, 5.0, 6.0}));
}

TEST_CASE("TensorTest-testPartialInvalidIndicesLength") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE_THROWS_AS(tensor.partial({0}, {1}), invalid_argument);
}

TEST_CASE("TensorTest-testPartialInvalidIndicesRange") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE_THROWS_AS(tensor.partial({0, 0}, {3, 1}), out_of_range);
}

TEST_CASE("TensorTest-testToStringMethod") {
    const Tensor tensor = Tensor({{1.0, 2.0}, {3.0, 4.0}});
    REQUIRE("Tensor(shape=[2, 2], data=[[1.0, 2.0], [3.0, 4.0]])" == tensor.to_string());
}
