//
// Created by Yiğit Demirşan on 21.02.2025.
//

#include "Tensor.h"
#include <stdexcept>
#include <numeric>
#include <sstream>
#include <functional>

Tensor::Tensor(const vector<vector<vector<float>>> &nested_data) {
    shape = infer_shape(nested_data);
    data = flatten(nested_data);
    strides = compute_strides(shape);
}

Tensor::Tensor(const vector<float> &flat_data, const vector<int> &shape) : data(flat_data), shape(shape) {
    if (compute_num_elements(shape) != flat_data.size()) {
        throw invalid_argument("Shape does not match the number of elements in data.");
    }
    strides = compute_strides(shape);
}

float Tensor::get(const vector<int> &indices) const {
    int flat_index = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        flat_index += indices[i] * strides[i];
    }
    return data[flat_index];
}

void Tensor::set(const vector<int> &indices, float value) {
    int flat_index = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        flat_index += indices[i] * strides[i];
    }
    data[flat_index] = value;
}

Tensor Tensor::reshape(const vector<int> &new_shape) const {
    if (compute_num_elements(new_shape) != compute_num_elements(shape)) {
        throw invalid_argument("Total number of elements must remain the same.");
    }
    return Tensor(data, new_shape);
}

Tensor Tensor::transpose(const vector<int> &axes) const {
    vector<int> new_shape(shape.size());
    for (size_t i = 0; i < shape.size(); ++i) {
        new_shape[i] = shape[axes[i]];
    }

    vector<float> new_data(compute_num_elements(new_shape));

    for (int i = 0; i < compute_num_elements(new_shape); ++i) {
        vector<int> new_indices = unflatten_index(i, compute_strides(new_shape));
        vector<int> original_indices(shape.size());
        for (size_t j = 0; j < shape.size(); ++j) {
            original_indices[axes[j]] = new_indices[j];
        }
        new_data[i] = get(original_indices);
    }

    return Tensor(new_data, new_shape);
}

Tensor Tensor::operator+(const Tensor &other) const {
    return elementwise_op(other, [](float a, float b) { return a + b; });
}

Tensor Tensor::operator-(const Tensor &other) const {
    return elementwise_op(other, [](float a, float b) { return a - b; });
}

Tensor Tensor::operator*(const Tensor &other) const {
    return elementwise_op(other, [](float a, float b) { return a * b; });
}

Tensor Tensor::dot(const Tensor &other) const {
    if (shape.back() != other.shape.front()) {
        throw invalid_argument("Shapes are not aligned for dot product.");
    }

    vector<int> result_shape = shape;
    result_shape.back() = other.shape.back();
    vector<float> result_data(compute_num_elements(result_shape), 0);

    for (int i = 0; i < compute_num_elements(result_shape); ++i) {
        vector<int> result_indices = unflatten_index(i, compute_strides(result_shape));
        float dot_product = 0;

        for (int k = 0; k < shape.back(); ++k) {
            vector<int> a_indices = result_indices;
            a_indices.back() = k;
            vector<int> b_indices = {k, result_indices.back()};

            dot_product += get(a_indices) * other.get(b_indices);
        }

        result_data[i] = dot_product;
    }

    return Tensor(result_data, result_shape);
}

string Tensor::to_string() const {
    ostringstream oss;
    oss << "Tensor(shape=[";
    for (size_t i = 0; i < shape.size(); ++i) {
        oss << shape[i] << (i + 1 < shape.size() ? ", " : "");
    }
    oss << "], data=[";
    for (size_t i = 0; i < data.size(); ++i) {
        oss << data[i] << (i + 1 < data.size() ? ", " : "");
    }
    oss << "])";
    return oss.str();
}

template <typename T>
vector<int> Tensor::infer_shape(const vector<T> &nested_data) {
    if constexpr (is_same<T, float>::value) {
        return {};
    } else {
        if (nested_data.empty()) {
            return {0};
        }
        vector<int> inner_shape = infer_shape(nested_data[0]);
        inner_shape.insert(inner_shape.begin(), nested_data.size());
        return inner_shape;
    }
}

template <typename T>
vector<float> Tensor::flatten(const vector<T> &nested_data) {
    if constexpr (is_same<T, float>::value) {
        return {nested_data};
    } else {
        vector<float> result;
        for (const auto &sublist : nested_data) {
            vector<float> subflattened = flatten(sublist);
            result.insert(result.end(), subflattened.begin(), subflattened.end());
        }
        return result;
    }
}

vector<int> Tensor::compute_strides(const vector<int> &shape) const {
    vector<int> strides(shape.size());
    int product = 1;
    for (int i = shape.size() - 1; i >= 0; --i) {
        strides[i] = product;
        product *= shape[i];
    }
    return strides;
}

int Tensor::compute_num_elements(const vector<int> &shape) const {
    return accumulate(shape.begin(), shape.end(), 1, multiplies<int>());
}

vector<int> Tensor::unflatten_index(int flat_index, const vector<int> &strides) const {
    vector<int> indices(strides.size());
    for (size_t i = 0; i < strides.size(); ++i) {
        indices[i] = flat_index / strides[i];
        flat_index %= strides[i];
    }
    return indices;
}

template <typename Op>
Tensor Tensor::elementwise_op(const Tensor &other, Op op) const {
    if (shape != other.shape) {
        throw invalid_argument("Shapes must match for element-wise operations.");
    }
    vector<float> result_data(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        result_data[i] = op(data[i], other.data[i]);
    }
    return Tensor(result_data, shape);
}
