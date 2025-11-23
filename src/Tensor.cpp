//
// Created by Yiğit Demirşan on 21.02.2025.
//

#include "Tensor.h"
#include <stdexcept>
#include <numeric>
#include <sstream>
#include <algorithm>

// --- Constructors ---

/**
 * Initializes the tensor with given data and shape.
 *
 * @param flatData  Nested lists representing the tensor data.
 * @param shape The shape of the tensor. If null, the shape is inferred from the data.
 */
Tensor::Tensor(const std::vector<float>& flatData, const std::vector<int>& shape)
        : data(flatData), shape(shape), strides(computeStrides(shape)) {
    if (computeNumberOfElements(shape) != static_cast<int>(flatData.size()))
        throw std::invalid_argument("Shape does not match data size.");
}

Tensor::Tensor(const std::vector<float>& nested_data)
        : shape({ static_cast<int>(nested_data.size()) }) {
    data = nested_data;
    strides = computeStrides(shape);
}

Tensor::Tensor(const std::vector<std::vector<float>>& nested_data)
        : shape(inferShape(nested_data)) {
    strides = computeStrides(shape);
}

Tensor::Tensor(const std::vector<std::vector<std::vector<float>>>& nested_data)
        : shape(inferShape(nested_data)) {
    strides = computeStrides(shape);
}

// --- Shape ---

std::vector<int> Tensor::getShape() const { return shape; }

/**
 * Computes the strides for each dimension based on the shape.
 *
 * @param shape Array representing the tensor shape.
 * @return Array representing the strides.
 */
std::vector<int> Tensor::computeStrides(const std::vector<int>& shape) {
    std::vector<int> strides(shape.size());
    int prod = 1;
    for (int i = static_cast<int>(shape.size()) - 1; i >= 0; --i) {
        strides[i] = prod;
        prod *= shape[i];
    }
    return strides;
}

/**
 * Computes the total number of elements in the tensor based on its shape.
 *
 * @param shape Array representing the tensor shape.
 * @return Total number of elements.
 */
int Tensor::computeNumberOfElements(const std::vector<int>& shape) {
    if (shape.empty()) return 0;
    return std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<int>());
}

/**
 * Infers the shape of the tensor from nested lists.
 *
 * @param data Nested lists representing the tensor data.
 * @return Array representing the shape.
 */
std::vector<int> Tensor::inferShape(const std::vector<float>& data) {
    return { static_cast<int>(data.size()) };
}
std::vector<int> Tensor::inferShape(const std::vector<std::vector<float>>& data) {
    if (data.empty()) return { 0 };
    return { static_cast<int>(data.size()), static_cast<int>(data[0].size()) };
}
std::vector<int> Tensor::inferShape(const std::vector<std::vector<std::vector<float>>>& data) {
    if (data.empty() || data[0].empty()) return { 0, 0, 0 };
    return { static_cast<int>(data.size()), static_cast<int>(data[0].size()), static_cast<int>(data[0][0].size()) };
}

/**
 * Validates that indices are within the valid range for each dimension.
 * @param indices Array of indices specifying the position.
 */
void Tensor::validateIndices(const std::vector<int>& indices) const {
    if (indices.size() != shape.size())
        throw std::out_of_range("Number of indices does not match number of dimensions.");
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] < 0 || indices[i] >= shape[i])
            throw std::out_of_range("Index out of bounds.");
    }
}

/**
 * Retrieves the value at the given indices.
 *
 * @param indices Array of indices specifying the position.
 * @return Value at the specified position.
 */
float Tensor::getValue(const std::vector<int>& indices) const {
    validateIndices(indices);
    int flat = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        flat += indices[i] * strides[i];
    }
    return data[flat];
}

/**
 * Sets the value at the given indices.
 *
 * @param indices Array of indices specifying the position.
 * @param value   Value to set at the specified position.
 */
void Tensor::setValue(const std::vector<int>& indices, float value) {
    validateIndices(indices);
    int flat = 0;
    for (size_t i = 0; i < indices.size(); ++i) flat += indices[i] * strides[i];
    data[flat] = value;
}

/**
 * Converts a flat index to multidimensional indices based on strides.
 *
 * @param flatIndex The flat index to convert.
 * @param shape   Array representing the strides.
 * @return Array of multi-dimensional indices.
 */
std::vector<int> Tensor::unflattenIndex(int flatIndex, const std::vector<int>& shape) {
    std::vector<int> strides = computeStrides(shape);
    std::vector<int> indices(shape.size());
    for (size_t i = 0; i < shape.size(); ++i) {
        indices[i] = flatIndex / strides[i];
        flatIndex %= strides[i];
    }
    return indices;
}

/**
 * Reshapes the tensor to the specified new shape.
 *
 * @param newShape Array representing the new shape.
 * @return New tensor with the specified shape.
 */
Tensor Tensor::reshape(const std::vector<int>& newShape) const {
    if (computeNumberOfElements(newShape) != computeNumberOfElements(shape))
        throw std::invalid_argument("Total number of elements must remain the same.");
    return Tensor{data, newShape};
}

/**
 * Transposes the tensor according to the specified axes.
 *
 * @param axes Array representing the order of axes. If null, reverses the axes.
 * @return New tensor with transposed axes.
 */
Tensor Tensor::transpose(const std::vector<int>& axes) const {
    std::vector<int> axes_order = axes;
    if (axes_order.empty()) {
        for (int i = static_cast<int>(shape.size()) - 1; i >= 0; --i) {
            axes_order.push_back(i);
        }
    }
    if (axes_order.size() != shape.size())
        throw std::invalid_argument("Invalid axes size for transpose.");
    std::vector<int> new_shape(shape.size());
    for (size_t i = 0; i < axes_order.size(); ++i) {
        new_shape[i] = shape[axes_order[i]];
    }
    std::vector<int> new_strides = computeStrides(new_shape);
    std::vector<float> new_data(data.size());
    for (int i = 0; i < computeNumberOfElements(new_shape); ++i) {
        std::vector<int> new_idx = unflattenIndex(i, new_shape);
        std::vector<int> orig_idx(shape.size());
        for (size_t j = 0; j < shape.size(); ++j) {
            orig_idx[axes_order[j]] = new_idx[j];
        }
        new_data[i] = getValue(orig_idx);
    }
    return Tensor{new_data, new_shape};
}

/**
 * Concatenates two tensors into a one.
 *
 * @param tensor 2nd tensor for concatenation.
 * @param dimension to concatenate.
 * @return Concatenated {@link Tensor}.
 */
Tensor Tensor::concat(const Tensor &tensor, int dimension) const {
    int startIndex = 1;
    int endIndex1 = 1;
    int endIndex2 = 1;
    for (int i = 0; i < this->shape.size(); i++) {
        if (i >= dimension) {
            endIndex1 *= this->shape[i];
            endIndex2 *= tensor.shape[i];
        } else {
            startIndex *= this->shape[i];
        }
    }
    std::vector<int> newShape(this->shape.size());
    for (int i = 0; i < this->shape.size(); i++) {
        if (i == dimension) {
            newShape[i] = this->shape[i] + tensor.shape[i];
        } else {
            newShape[i] = this->shape[i];
        }
    }
    std::vector<float> newList(startIndex * (endIndex1 + endIndex2));
    int k = 0;
    for (int i = 0; i < startIndex; i++) {
        for (int j = 0; j < endIndex1; j++) {
            newList[k++] = this->data[i * endIndex1 + j];
        }
        for (int j = 0; j < endIndex2; j++) {
            newList[k++] = tensor.data[i * endIndex2 + j];
        }
    }
    return Tensor{newList, newShape};
}

/**
 * Returns the sub-{@link Tensor} taking the given dimensions.
 *
 * @return a sub-{@link Tensor}.
 */
Tensor Tensor::get(const std::vector<int> &dimensions) const {
    std::vector<int> newShape(this->shape.size() - dimensions.size());
    newShape.insert(newShape.begin(), this->shape.begin() + static_cast<int>(dimensions.size()), this->shape.end());
    int i = 0, start = 0, end = static_cast<int>(data.size());
    do {
        int parts = (end - start) / this->shape[i];
        start += parts * dimensions[i];
        end = start + parts;
        i++;
    } while (i < dimensions.size());
    std::vector<float> sublist;
    sublist.insert(sublist.begin(), this->data.begin() + start, this->data.begin() + end);
    return Tensor{sublist,  newShape};
}

/**
 * Determines the broadcasted shape of two tensors.
 *
 * @param shape1 Array representing the first tensor shape.
 * @param shape2 Array representing the second tensor shape.
 * @return Array representing the broadcasted shape.
 */
std::vector<int> Tensor::broadcastShape(const std::vector<int>& shape1, const std::vector<int>& shape2) {
    size_t ndim = std::max(shape1.size(), shape2.size());
    std::vector<int> s1 = shape1, s2 = shape2;
    s1.insert(s1.begin(), ndim - s1.size(), 1);
    s2.insert(s2.begin(), ndim - s2.size(), 1);
    std::vector<int> out(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        if (s1[i] == s2[i] || s2[i] == 1) out[i] = s1[i];
        else if (s1[i] == 1) out[i] = s2[i];
        else throw std::invalid_argument("Shapes are not broadcastable");
    }
    return out;
}

float Tensor::getBroadcasted(const std::vector<int>& indices, const std::vector<int>& expanded_shape) const {
    std::vector<int> offset_indices;
    int dim_offset = static_cast<int>(expanded_shape.size() - shape.size());
    for (size_t i = 0; i < shape.size(); ++i) {
        if (shape[i] == 1)
            offset_indices.push_back(0);
        else
            offset_indices.push_back(indices[dim_offset + i]);
    }
    return getValue(offset_indices);
}

/**
 * Broadcasts the tensor to the specified target shape.
 *
 * @param targetShape Array representing the target shape.
 * @return New tensor with the target shape.
 */
Tensor Tensor::broadcastTo(const std::vector<int>& targetShape) const {
    std::vector<int> expanded_shape = shape;
    expanded_shape.insert(expanded_shape.begin(), targetShape.size() - shape.size(), 1);
    for (size_t i = 0; i < targetShape.size(); ++i) {
        if (!(expanded_shape[i] == targetShape[i] || expanded_shape[i] == 1))
            throw std::invalid_argument("Cannot broadcast to target shape.");
    }
    int total = computeNumberOfElements(targetShape);
    std::vector<float> new_data;
    for (int idx = 0; idx < total; ++idx) {
        auto idxs = unflattenIndex(idx, targetShape);
        new_data.push_back(getBroadcasted(idxs, targetShape));
    }
    return Tensor{new_data, targetShape};
}

/**
 * Adds two tensors element-wise with broadcasting.
 *
 * @param other The other tensor to add.
 * @return New tensor with the result of the addition.
 */
Tensor Tensor::add(const Tensor& other) const {
    auto out_shape = broadcastShape(shape, other.shape);
    auto t1 = broadcastTo(out_shape), t2 = other.broadcastTo(out_shape);
    int total = computeNumberOfElements(out_shape);
    std::vector<float> result(total);
    for (int i = 0; i < total; ++i)
        result[i] = t1.data[i] + t2.data[i];
    return Tensor{result, out_shape};
}

/**
 * Subtracts one tensor from another element-wise with broadcasting.
 *
 * @param other The other tensor to subtract.
 * @return New tensor with the result of the subtraction.
 */
Tensor Tensor::subtract(const Tensor& other) const {
    auto out_shape = broadcastShape(shape, other.shape);
    auto t1 = broadcastTo(out_shape), t2 = other.broadcastTo(out_shape);
    int total = computeNumberOfElements(out_shape);
    std::vector<float> result(total);
    for (int i = 0; i < total; ++i)
        result[i] = t1.data[i] - t2.data[i];
    return Tensor{result, out_shape};
}

/**
 * Multiplies two tensors element-wise with broadcasting.
 *
 * @param other The other tensor to multiply.
 * @return New tensor with the result of the multiplication.
 */
Tensor Tensor::hadamardProduct(const Tensor& other) const {
    auto out_shape = broadcastShape(shape, other.shape);
    auto t1 = broadcastTo(out_shape), t2 = other.broadcastTo(out_shape);
    int total = computeNumberOfElements(out_shape);
    std::vector<float> result(total);
    for (int i = 0; i < total; ++i)
        result[i] = t1.data[i] * t2.data[i];
    return Tensor{result, out_shape};
}

/**
 * Performs matrix multiplication (batched if necessary).
 * For tensors of shape (..., M, K) and (..., K, N), returns (..., M, N).
 *
 * @param other Tensor with shape compatible for matrix multiplication.
 * @return Tensor resulting from matrix multiplication.
 */
Tensor Tensor::multiply(const Tensor& other) const {
    // Both must be at least 2D (batch, m, k) and (batch, k, n)
    if (shape.size() < 2 || other.shape.size() < 2)
        throw std::invalid_argument("matmul only supports tensors with at least 2 dimensions.");
    int m = shape[shape.size() - 2];
    int k1 = shape[shape.size() - 1];
    int k2 = other.shape[other.shape.size() - 2];
    int n = other.shape[other.shape.size() - 1];
    if (k1 != k2)
        throw std::invalid_argument("matmul: inner dimensions do not match.");
    std::vector<int> batch_shape1(shape.begin(), shape.end() - 2);
    std::vector<int> batch_shape2(other.shape.begin(), other.shape.end() - 2);
    std::vector<int> batch_shape = broadcastShape(batch_shape1, batch_shape2);
    std::vector<int> shapeA = batch_shape;
    shapeA.push_back(m);
    shapeA.push_back(k1);
    std::vector<int> shapeB = batch_shape;
    shapeB.push_back(k2);
    shapeB.push_back(n);
    Tensor A = broadcastTo(shapeA);
    Tensor B = other.broadcastTo(shapeB);
    std::vector<int> result_shape = batch_shape;
    result_shape.push_back(m);
    result_shape.push_back(n);
    int result_total = computeNumberOfElements(result_shape);
    std::vector<float> result_data(result_total);
    for (int idx = 0; idx < result_total; ++idx) {
        std::vector<int> indices = unflattenIndex(idx, result_shape);
        std::vector<int> batch_indices(indices.begin(), indices.end() - 2);
        int row = indices[indices.size() - 2];
        int col = indices[indices.size() - 1];
        float sum = 0.0f;
        for (int kk = 0; kk < k1; ++kk) {
            std::vector<int> a_idx = batch_indices;
            a_idx.push_back(row); a_idx.push_back(kk);
            std::vector<int> b_idx = batch_indices;
            b_idx.push_back(kk); b_idx.push_back(col);
            sum += A.getValue(a_idx) * B.getValue(b_idx);
        }
        result_data[idx] = sum;
    }
    return Tensor{result_data, result_shape};
}

/**
 * Extracts a sub-tensor from the given start indices to the end indices.
 *
 * @param startIndices Array specifying the start indices for each dimension.
 * @param endIndices   Array specifying the end indices (exclusive) for each dimension.
 * @return A new Tensor containing the extracted sub-tensor.
 */
Tensor Tensor::partial(const std::vector<int>& startIndices, const std::vector<int>& endIndices) const {
    if (startIndices.size() != shape.size() || endIndices.size() != shape.size())
        throw std::invalid_argument("partial: indices must match number of dimensions.");
    std::vector<int> new_shape(shape.size());
    for (size_t i = 0; i < shape.size(); ++i) {
        if (endIndices[i] < startIndices[i]) throw std::invalid_argument("partial: end < start");
        new_shape[i] = endIndices[i] - startIndices[i];
    }
    int total = computeNumberOfElements(new_shape);
    std::vector<float> new_data(total);
    for (int idx = 0; idx < total; ++idx) {
        std::vector<int> sub_idx = unflattenIndex(idx, new_shape);
        std::vector<int> orig_idx(shape.size());
        for (size_t j = 0; j < shape.size(); ++j) {
            orig_idx[j] = startIndices[j] + sub_idx[j];
        }
        new_data[idx] = getValue(orig_idx);
    }
    return Tensor{new_data, new_shape};
}

/**
 * Returns a string representation of the tensor.
 *
 * @return String representing the tensor.
 */
std::string Tensor::to_string() const {
    std::ostringstream oss;
    oss << "Tensor(shape=[";
    for (size_t i = 0; i < shape.size(); ++i) {
        oss << shape[i];
        if (i + 1 < shape.size()) oss << ", ";
    }
    oss << "], data=[";
    for (size_t i = 0; i < data.size(); ++i) {
        oss << data[i];
        if (i + 1 < data.size()) oss << ", ";
    }
    oss << "])";
    return oss.str();
}
