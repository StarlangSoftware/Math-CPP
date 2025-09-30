//
// Created by Yiğit Demirşan on 21.02.2025.
//

#include "Tensor.h"
#include <stdexcept>
#include <numeric>
#include <sstream>
#include <algorithm>

// --- Constructors ---

Tensor::Tensor(const std::vector<float>& flat_data, const std::vector<int>& shape)
        : data_(flat_data), shape_(shape), strides_(compute_strides(shape)) {
    if (compute_num_elements(shape) != (int)flat_data.size())
        throw std::invalid_argument("Shape does not match data size.");
}

Tensor::Tensor(const std::vector<float>& nested_data)
        : shape_({ (int)nested_data.size() }) {
    data_ = nested_data;
    strides_ = compute_strides(shape_);
}

Tensor::Tensor(const std::vector<std::vector<float>>& nested_data)
        : shape_(infer_shape(nested_data)) {
    flatten_helper(nested_data, data_);
    strides_ = compute_strides(shape_);
}

Tensor::Tensor(const std::vector<std::vector<std::vector<float>>>& nested_data)
        : shape_(infer_shape(nested_data)) {
    flatten_helper(nested_data, data_);
    strides_ = compute_strides(shape_);
}

// --- Shape and strides ---

std::vector<int> Tensor::shape() const { return shape_; }
std::vector<int> Tensor::strides() const { return strides_; }

// --- Internal helpers ---

std::vector<int> Tensor::compute_strides(const std::vector<int>& shape) {
    std::vector<int> strides(shape.size());
    int prod = 1;
    for (int i = (int)shape.size() - 1; i >= 0; --i) {
        strides[i] = prod;
        prod *= shape[i];
    }
    return strides;
}

int Tensor::compute_num_elements(const std::vector<int>& shape) {
    if (shape.empty()) return 0;
    return std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<int>());
}

// Scalar version: only for float (or things convertible to float)
template<typename T>
typename std::enable_if<std::is_arithmetic<T>::value, void>::type
flatten_helper(const T& input, std::vector<float>& out) {
    out.push_back(static_cast<float>(input));
}

// Vector version: for any vector type
template<typename T>
typename std::enable_if<!std::is_arithmetic<T>::value, void>::type
flatten_helper(const T& input, std::vector<float>& out) {
    for (const auto& element : input) {
        flatten_helper(element, out);
    }
}

std::vector<int> Tensor::infer_shape(const std::vector<float>& data) {
    return { (int)data.size() };
}
std::vector<int> Tensor::infer_shape(const std::vector<std::vector<float>>& data) {
    if (data.empty()) return { 0 };
    return { (int)data.size(), (int)data[0].size() };
}
std::vector<int> Tensor::infer_shape(const std::vector<std::vector<std::vector<float>>>& data) {
    if (data.empty() || data[0].empty()) return { 0, 0, 0 };
    return { (int)data.size(), (int)data[0].size(), (int)data[0][0].size() };
}

// --- Indexing and validation ---

void Tensor::validate_indices(const std::vector<int>& indices) const {
    if (indices.size() != shape_.size())
        throw std::out_of_range("Number of indices does not match number of dimensions.");
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] < 0 || indices[i] >= shape_[i])
            throw std::out_of_range("Index out of bounds.");
    }
}

float Tensor::get(const std::vector<int>& indices) const {
    validate_indices(indices);
    int flat = 0;
    for (size_t i = 0; i < indices.size(); ++i) flat += indices[i] * strides_[i];
    return data_[flat];
}
void Tensor::set(const std::vector<int>& indices, float value) {
    validate_indices(indices);
    int flat = 0;
    for (size_t i = 0; i < indices.size(); ++i) flat += indices[i] * strides_[i];
    data_[flat] = value;
}

std::vector<int> Tensor::unflatten_index(int flat_index, const std::vector<int>& shape) {
    std::vector<int> strides = compute_strides(shape);
    std::vector<int> indices(shape.size());
    for (size_t i = 0; i < shape.size(); ++i) {
        indices[i] = flat_index / strides[i];
        flat_index %= strides[i];
    }
    return indices;
}

// --- Reshape and transpose ---

Tensor Tensor::reshape(const std::vector<int>& new_shape) const {
    if (compute_num_elements(new_shape) != compute_num_elements(shape_))
        throw std::invalid_argument("Total number of elements must remain the same.");
    return Tensor(data_, new_shape);
}

Tensor Tensor::transpose(const std::vector<int>& axes) const {
    std::vector<int> axes_order = axes;
    if (axes_order.empty()) {
        for (int i = (int)shape_.size() - 1; i >= 0; --i) axes_order.push_back(i);
    }
    if (axes_order.size() != shape_.size())
        throw std::invalid_argument("Invalid axes size for transpose.");
    std::vector<int> new_shape(shape_.size());
    for (size_t i = 0; i < axes_order.size(); ++i) new_shape[i] = shape_[axes_order[i]];
    std::vector<int> new_strides = compute_strides(new_shape);
    std::vector<float> new_data(data_.size());
    for (int i = 0; i < compute_num_elements(new_shape); ++i) {
        std::vector<int> new_idx = unflatten_index(i, new_shape);
        std::vector<int> orig_idx(shape_.size());
        for (size_t j = 0; j < shape_.size(); ++j) orig_idx[axes_order[j]] = new_idx[j];
        new_data[i] = get(orig_idx);
    }
    return Tensor(new_data, new_shape);
}

// --- Broadcasting ---

std::vector<int> Tensor::broadcast_shape(const std::vector<int>& shape1, const std::vector<int>& shape2) {
    size_t ndim = std::max(shape1.size(), shape2.size());
    std::vector<int> s1 = shape1, s2 = shape2;
    s1.insert(s1.begin(), ndim - s1.size(), 1);
    s2.insert(s2.begin(), ndim - s2.size(), 1);
    std::vector<int> out(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        if (s1[i] == s2[i]) out[i] = s1[i];
        else if (s1[i] == 1) out[i] = s2[i];
        else if (s2[i] == 1) out[i] = s1[i];
        else throw std::invalid_argument("Shapes are not broadcastable");
    }
    return out;
}

float Tensor::get_broadcasted(const std::vector<int>& indices, const std::vector<int>& expanded_shape) const {
    std::vector<int> offset_indices;
    int dim_offset = expanded_shape.size() - shape_.size();
    for (size_t i = 0; i < shape_.size(); ++i) {
        if (shape_[i] == 1)
            offset_indices.push_back(0);
        else
            offset_indices.push_back(indices[dim_offset + i]);
    }
    return get(offset_indices);
}

Tensor Tensor::broadcast_to(const std::vector<int>& target_shape) const {
    std::vector<int> expanded_shape = shape_;
    expanded_shape.insert(expanded_shape.begin(), target_shape.size() - shape_.size(), 1);
    for (size_t i = 0; i < target_shape.size(); ++i) {
        if (!(expanded_shape[i] == target_shape[i] || expanded_shape[i] == 1))
            throw std::invalid_argument("Cannot broadcast to target shape.");
    }
    int total = compute_num_elements(target_shape);
    std::vector<float> new_data;
    for (int idx = 0; idx < total; ++idx) {
        auto idxs = unflatten_index(idx, target_shape);
        new_data.push_back(get_broadcasted(idxs, target_shape));
    }
    return Tensor(new_data, target_shape);
}

// --- Elementwise ops ---

Tensor Tensor::add(const Tensor& other) const {
    auto out_shape = broadcast_shape(shape_, other.shape_);
    auto t1 = broadcast_to(out_shape), t2 = other.broadcast_to(out_shape);
    std::vector<float> result;
    int total = compute_num_elements(out_shape);
    for (int i = 0; i < total; ++i)
        result.push_back(t1.data_[i] + t2.data_[i]);
    return Tensor(result, out_shape);
}
Tensor Tensor::subtract(const Tensor& other) const {
    auto out_shape = broadcast_shape(shape_, other.shape_);
    auto t1 = broadcast_to(out_shape), t2 = other.broadcast_to(out_shape);
    std::vector<float> result;
    int total = compute_num_elements(out_shape);
    for (int i = 0; i < total; ++i)
        result.push_back(t1.data_[i] - t2.data_[i]);
    return Tensor(result, out_shape);
}
Tensor Tensor::hadamard_product(const Tensor& other) const {
    auto out_shape = broadcast_shape(shape_, other.shape_);
    auto t1 = broadcast_to(out_shape), t2 = other.broadcast_to(out_shape);
    std::vector<float> result;
    int total = compute_num_elements(out_shape);
    for (int i = 0; i < total; ++i)
        result.push_back(t1.data_[i] * t2.data_[i]);
    return Tensor(result, out_shape);
}

// --- Matrix multiplication (batched matmul) ---

Tensor Tensor::matmul(const Tensor& other) const {
    // Both must be at least 2D (batch, m, k) and (batch, k, n)
    if (shape_.size() < 2 || other.shape_.size() < 2)
        throw std::invalid_argument("matmul only supports tensors with at least 2 dimensions.");
    int m = shape_[shape_.size() - 2];
    int k1 = shape_[shape_.size() - 1];
    int k2 = other.shape_[other.shape_.size() - 2];
    int n = other.shape_[other.shape_.size() - 1];
    if (k1 != k2)
        throw std::invalid_argument("matmul: inner dimensions do not match.");

    std::vector<int> batch_shape1(shape_.begin(), shape_.end() - 2);
    std::vector<int> batch_shape2(other.shape_.begin(), other.shape_.end() - 2);
    std::vector<int> batch_shape = broadcast_shape(batch_shape1, batch_shape2);

    std::vector<int> shapeA = batch_shape;
    shapeA.push_back(m);
    shapeA.push_back(k1);

    std::vector<int> shapeB = batch_shape;
    shapeB.push_back(k2);
    shapeB.push_back(n);

    Tensor A = broadcast_to(shapeA);
    Tensor B = other.broadcast_to(shapeB);


    std::vector<int> result_shape = batch_shape;
    result_shape.push_back(m);
    result_shape.push_back(n);

    int result_total = compute_num_elements(result_shape);
    std::vector<float> result_data(result_total);

    for (int idx = 0; idx < result_total; ++idx) {
        std::vector<int> indices = unflatten_index(idx, result_shape);
        std::vector<int> batch_indices(indices.begin(), indices.end() - 2);
        int row = indices[indices.size() - 2];
        int col = indices[indices.size() - 1];
        float sum = 0.0f;
        for (int kk = 0; kk < k1; ++kk) {
            std::vector<int> a_idx = batch_indices;
            a_idx.push_back(row); a_idx.push_back(kk);
            std::vector<int> b_idx = batch_indices;
            b_idx.push_back(kk); b_idx.push_back(col);
            sum += A.get(a_idx) * B.get(b_idx);
        }
        result_data[idx] = sum;
    }
    return Tensor(result_data, result_shape);
}

// --- Partial/slice ---

Tensor Tensor::partial(const std::vector<int>& start_indices, const std::vector<int>& end_indices) const {
    if (start_indices.size() != shape_.size() || end_indices.size() != shape_.size())
        throw std::invalid_argument("partial: indices must match number of dimensions.");
    std::vector<int> new_shape(shape_.size());
    for (size_t i = 0; i < shape_.size(); ++i) {
        if (end_indices[i] < start_indices[i]) throw std::invalid_argument("partial: end < start");
        new_shape[i] = end_indices[i] - start_indices[i];
    }
    int total = compute_num_elements(new_shape);
    std::vector<float> new_data(total);
    for (int idx = 0; idx < total; ++idx) {
        std::vector<int> sub_idx = unflatten_index(idx, new_shape);
        std::vector<int> orig_idx(shape_.size());
        for (size_t j = 0; j < shape_.size(); ++j)
            orig_idx[j] = start_indices[j] + sub_idx[j];
        new_data[idx] = get(orig_idx);
    }
    return Tensor(new_data, new_shape);
}

// --- String representation ---

std::string Tensor::to_string() const {
    std::ostringstream oss;
    oss << "Tensor(shape=[";
    for (size_t i = 0; i < shape_.size(); ++i) {
        oss << shape_[i];
        if (i + 1 < shape_.size()) oss << ", ";
    }
    oss << "], data=[";
    for (size_t i = 0; i < data_.size(); ++i) {
        oss << data_[i];
        if (i + 1 < data_.size()) oss << ", ";
    }
    oss << "])";
    return oss.str();
}
