//
// Created by Yiğit Demirşan on 21.02.2025.
//

#ifndef MATH_TENSOR_H
#define MATH_TENSOR_H

#include <string>
#include <vector>
#include <string>
#include <tuple>
#include <initializer_list>

class Tensor {
public:
    // Constructors
    Tensor(const std::vector<float>& flat_data, const std::vector<int>& shape);
    Tensor(const std::vector<std::vector<std::vector<float>>>& nested_data); // For 3D data
    Tensor(const std::vector<std::vector<float>>& nested_data);              // For 2D data
    Tensor(const std::vector<float>& nested_data);                           // For 1D data

    // Element access
    float get_value(const std::vector<int>& indices) const;
    void set(const std::vector<int>& indices, float value);

    Tensor get(const std::vector<int>& indices) const;
    Tensor concat(Tensor& tensor, int dimension);

    // Shape and strides
    std::vector<int> shape() const;
    std::vector<int> strides() const;

    // Reshape and transpose
    Tensor reshape(const std::vector<int>& new_shape) const;
    Tensor transpose(const std::vector<int>& axes = {}) const;

    // Broadcasting and elementwise ops
    Tensor broadcast_to(const std::vector<int>& target_shape) const;
    Tensor add(const Tensor& other) const;
    Tensor subtract(const Tensor& other) const;
    Tensor hadamard_product(const Tensor& other) const;

    // Matrix multiplication
    Tensor matmul(const Tensor& other) const;

    // Slicing
    Tensor partial(const std::vector<int>& start_indices, const std::vector<int>& end_indices) const;

    // Representation
    std::string to_string() const;

private:
    std::vector<float> data_;
    std::vector<int> shape_;
    std::vector<int> strides_;

    // Utilities
    static std::vector<int> infer_shape(const std::vector<float>& data);
    static std::vector<int> infer_shape(const std::vector<std::vector<float>>& data);
    static std::vector<int> infer_shape(const std::vector<std::vector<std::vector<float>>>& data);

    template <typename T>
    static void flatten_helper(const T& input, std::vector<float>& out);

    static std::vector<int> compute_strides(const std::vector<int>& shape);
    static int compute_num_elements(const std::vector<int>& shape);
    void validate_indices(const std::vector<int>& indices) const;
    static std::vector<int> unflatten_index(int flat_index, const std::vector<int>& shape);

    // Broadcasting
    static std::vector<int> broadcast_shape(const std::vector<int>& shape1, const std::vector<int>& shape2);

    // Internal for broadcast_to
    float get_broadcasted(const std::vector<int>& indices, const std::vector<int>& expanded_shape) const;
};

#endif //MATH_TENSOR_H
