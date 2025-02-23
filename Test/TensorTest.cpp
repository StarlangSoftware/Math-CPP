//
// Created by Yiğit Demirşan on 23.02.2025.
//

#include "src/Tensor.h"
#include <gtest/gtest.h>

class TensorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialization before each test if needed
    }

    void TearDown() override {
        // Cleanup after each test if needed
    }
};

// Test tensor initialization with nested data
TEST_F(TensorTest, InitializationFromNestedData) {
vector<vector<vector<float>>> nested_data = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};
Tensor tensor(nested_data);

EXPECT_EQ(tensor.get({0, 0, 0}), 1);
EXPECT_EQ(tensor.get({1, 1, 1}), 8);
}

// Test tensor initialization from flat data
TEST_F(TensorTest, InitializationFromFlatData) {
vector<float> flat_data = {1, 2, 3, 4, 5, 6};
vector<int> shape = {2, 3};
Tensor tensor(flat_data, shape);

EXPECT_EQ(tensor.get({0, 0}), 1);
EXPECT_EQ(tensor.get({1, 2}), 6);
}

// Test get and set methods
TEST_F(TensorTest, GetSet) {
vector<float> flat_data = {1, 2, 3, 4, 5, 6};
vector<int> shape = {2, 3};
Tensor tensor(flat_data, shape);

tensor.set({0, 1}, 10);
EXPECT_EQ(tensor.get({0, 1}), 10);
}

// Test reshape method
TEST_F(TensorTest, Reshape) {
vector<float> flat_data = {1, 2, 3, 4};
vector<int> shape = {2, 2};
Tensor tensor(flat_data, shape);

Tensor reshaped = tensor.reshape({4});
EXPECT_EQ(reshaped.get({2}), 3);
}

// Test transpose method
TEST_F(TensorTest, Transpose) {
vector<float> flat_data = {1, 2, 3, 4, 5, 6};
vector<int> shape = {2, 3};
Tensor tensor(flat_data, shape);

Tensor transposed = tensor.transpose({1, 0});
EXPECT_EQ(transposed.get({0, 1}), 3);
EXPECT_EQ(transposed.get({2, 1}), 6);
}

// Test addition operator
TEST_F(TensorTest, Addition) {
vector<float> flat_data1 = {1, 2, 3, 4};
vector<float> flat_data2 = {5, 6, 7, 8};
vector<int> shape = {2, 2};

Tensor tensor1(flat_data1, shape);
Tensor tensor2(flat_data2, shape);

Tensor result = tensor1 + tensor2;

EXPECT_EQ(result.get({0, 0}), 6);
EXPECT_EQ(result.get({1, 1}), 12);
}

// Test subtraction operator
TEST_F(TensorTest, Subtraction) {
vector<float> flat_data1 = {5, 6, 7, 8};
vector<float> flat_data2 = {1, 2, 3, 4};
vector<int> shape = {2, 2};

Tensor tensor1(flat_data1, shape);
Tensor tensor2(flat_data2, shape);

Tensor result = tensor1 - tensor2;

EXPECT_EQ(result.get({0, 0}), 4);
EXPECT_EQ(result.get({1, 1}), 4);
}

// Test multiplication operator
TEST_F(TensorTest, Multiplication) {
vector<float> flat_data1 = {1, 2, 3, 4};
vector<float> flat_data2 = {2, 3, 4, 5};
vector<int> shape = {2, 2};

Tensor tensor1(flat_data1, shape);
Tensor tensor2(flat_data2, shape);

Tensor result = tensor1 * tensor2;

EXPECT_EQ(result.get({0, 0}), 2);
EXPECT_EQ(result.get({1, 1}), 20);
}

// Test dot product
TEST_F(TensorTest, DotProduct) {
vector<float> flat_data1 = {1, 2, 3, 4};
vector<float> flat_data2 = {1, 2, 3, 4};
vector<int> shape1 = {2, 2};
vector<int> shape2 = {2, 2};

Tensor tensor1(flat_data1, shape1);
Tensor tensor2(flat_data2, shape2);

Tensor result = tensor1.dot(tensor2);

EXPECT_EQ(result.get({0, 0}), 7);
EXPECT_EQ(result.get({1, 1}), 25);
}

// Test to_string method
TEST_F(TensorTest, ToString) {
vector<float> flat_data = {1, 2, 3, 4};
vector<int> shape = {2, 2};
Tensor tensor(flat_data, shape);

string tensor_str = tensor.to_string();
EXPECT_TRUE(tensor_str.find("Tensor(shape=[2, 2], data=[1, 2, 3, 4])") != string::npos);
}

// Test invalid reshape
TEST_F(TensorTest, InvalidReshape) {
vector<float> flat_data = {1, 2, 3, 4};
vector<int> shape = {2, 2};
Tensor tensor(flat_data, shape);

EXPECT_THROW(tensor.reshape({3, 2}), invalid_argument);
}

// Test invalid dot product
TEST_F(TensorTest, InvalidDotProduct) {
vector<float> flat_data1 = {1, 2, 3, 4};
vector<float> flat_data2 = {1, 2, 3, 4};
vector<int> shape1 = {2, 2};
vector<int> shape2 = {3, 2};

Tensor tensor1(flat_data1, shape1);
Tensor tensor2(flat_data2, shape2);

EXPECT_THROW(tensor1.dot(tensor2), invalid_argument);
}

