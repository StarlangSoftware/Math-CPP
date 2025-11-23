//
// Created by Yiğit Demirşan on 21.02.2025.
//

#ifndef MATH_TENSOR_H
#define MATH_TENSOR_H

#include <string>
#include <vector>

class Tensor {
public:
    // Constructors
    Tensor(const std::vector<float>& flatData, const std::vector<int>& shape);
    explicit Tensor(const std::vector<std::vector<std::vector<float>>>& nested_data); // For 3D data
    explicit Tensor(const std::vector<std::vector<float>>& nested_data);              // For 2D data
    explicit Tensor(const std::vector<float>& nested_data);                           // For 1D data

    // Element access
    float getValue(const std::vector<int>& indices) const;
    void setValue(const std::vector<int>& indices, float value);

    // Shape and strides
    std::vector<int> getShape() const;

    // Reshape and transpose
    Tensor reshape(const std::vector<int>& newShape) const;
    Tensor transpose(const std::vector<int>& axes = {}) const;

    Tensor concat(const Tensor& other, int dimension) const;
    Tensor get(const std::vector<int>& dimensions) const;

    // Broadcasting and elementwise ops
    Tensor broadcastTo(const std::vector<int>& targetShape) const;
    Tensor add(const Tensor& other) const;
    Tensor subtract(const Tensor& other) const;
    Tensor hadamardProduct(const Tensor& other) const;

    // Matrix multiplication
    Tensor multiply(const Tensor& other) const;

    // Slicing
    Tensor partial(const std::vector<int>& startIndices, const std::vector<int>& endIndices) const;

    // Representation
    std::string to_string() const;

private:
    std::vector<float> data;
    std::vector<int> shape;
    std::vector<int> strides;

    // Utilities
    static std::vector<int> inferShape(const std::vector<float>& data);
    static std::vector<int> inferShape(const std::vector<std::vector<float>>& data);
    static std::vector<int> inferShape(const std::vector<std::vector<std::vector<float>>>& data);

    static std::vector<int> computeStrides(const std::vector<int>& shape);
    static int computeNumberOfElements(const std::vector<int>& shape);
    void validateIndices(const std::vector<int>& indices) const;
    static std::vector<int> unflattenIndex(int flatIndex, const std::vector<int>& shape);

    // Broadcasting
    static std::vector<int> broadcastShape(const std::vector<int>& shape1, const std::vector<int>& shape2);

    // Internal for broadcast_to
    float getBroadcasted(const std::vector<int>& indices, const std::vector<int>& expanded_shape) const;
};

#endif //MATH_TENSOR_H
