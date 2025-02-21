//
// Created by Yiğit Demirşan on 21.02.2025.
//

#ifndef MATH_TENSOR_H
#define MATH_TENSOR_H

#include <vector>

using namespace std;

class Tensor {
public:
    Tensor(const vector<vector<vector<float>>> &nested_data);
    Tensor(const vector<float> &flat_data, const vector<int> &shape);

    float get(const vector<int> &indices) const;
    void set(const vector<int> &indices, float value);

    Tensor reshape(const vector<int> &new_shape) const;
    Tensor transpose(const vector<int> &axes) const;

    Tensor operator+(const Tensor &other) const;
    Tensor operator-(const Tensor &other) const;
    Tensor operator*(const Tensor &other) const;
    Tensor dot(const Tensor &other) const;

    auto to_string() const;

private:
    vector<float> data;
    vector<int> shape;
    vector<int> strides;

    template <typename T>
    vector<int> infer_shape(const vector<T> &nested_data);

    template <typename T>
    vector<float> flatten(const vector<T> &nested_data);

    vector<int> compute_strides(const vector<int> &shape) const;
    int compute_num_elements(const vector<int> &shape) const;
    vector<int> unflatten_index(int flat_index, const vector<int> &strides) const;

    template <typename Op>
    Tensor elementwise_op(const Tensor &other, Op op) const;
};


#endif //MATH_TENSOR_H
