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
    [[nodiscard]] float getValue(const std::vector<int>& indices) const;
    void setValue(const std::vector<int>& indices, float value);

    // Shape and strides
    [[nodiscard]] std::vector<int> getShape() const;

    // Reshape and transpose
    [[nodiscard]] Tensor reshape(const std::vector<int>& newShape) const;
    [[nodiscard]] Tensor transpose(const std::vector<int>& axes = {}) const;

    [[nodiscard]] Tensor concat(const Tensor& other, int dimension) const;
    [[nodiscard]] Tensor get(const std::vector<int>& dimensions) const;

    // Broadcasting and elementwise ops
    [[nodiscard]] Tensor broadcastTo(const std::vector<int>& targetShape) const;
    [[nodiscard]] Tensor add(const Tensor& other) const;
    [[nodiscard]] Tensor subtract(const Tensor& other) const;
    [[nodiscard]] Tensor hadamardProduct(const Tensor& other) const;

    // Matrix multiplication
    [[nodiscard]] Tensor multiply(const Tensor& other) const;

    // Slicing
    [[nodiscard]] Tensor partial(const std::vector<int>& startIndices, const std::vector<int>& endIndices) const;

    // Representation
    [[nodiscard]] std::string to_string() const;

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
    [[nodiscard]] float getBroadcasted(const std::vector<int>& indices, const std::vector<int>& expanded_shape) const;
};

#endif //MATH_TENSOR_H
