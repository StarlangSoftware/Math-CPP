//
// Created by Yiğit Demirşan on 21.02.2025.
//

#ifndef MATH_TENSOR_H
#define MATH_TENSOR_H

#include <string>
#include <vector>

using namespace std;

class Tensor {
public:
    // Constructors
    Tensor(const vector<double>& flatData, const vector<int>& shape);
    explicit Tensor(const vector<vector<vector<double>>>& nested_data); // For 3D data
    explicit Tensor(const vector<vector<double>>& nested_data);              // For 2D data
    explicit Tensor(const vector<double>& nested_data);                           // For 1D data

    // Element access
    [[nodiscard]] double getValue(const vector<int>& indices) const;
    void setValue(const vector<int>& indices, double value);

    // Shape and strides
    [[nodiscard]] vector<int> getShape() const;

    [[nodiscard]] vector<double> getData() const;

    // Reshape and transpose
    [[nodiscard]] Tensor reshape(const vector<int>& newShape) const;
    [[nodiscard]] Tensor transpose(const vector<int>& axes = {}) const;

    [[nodiscard]] Tensor concat(const Tensor& other, int dimension) const;
    [[nodiscard]] Tensor get(const vector<int>& dimensions) const;

    // Broadcasting and elementwise ops
    [[nodiscard]] Tensor broadcastTo(const vector<int>& targetShape) const;
    [[nodiscard]] Tensor add(const Tensor& other) const;
    [[nodiscard]] Tensor subtract(const Tensor& other) const;
    [[nodiscard]] Tensor hadamardProduct(const Tensor& other) const;

    // Matrix multiplication
    [[nodiscard]] Tensor multiply(const Tensor& other) const;

    // Slicing
    [[nodiscard]] Tensor partial(const vector<int>& startIndices, const vector<int>& endIndices) const;

    // Representation
    [[nodiscard]] std::string to_string() const;

private:
    vector<double> data;
    vector<int> shape;
    vector<int> strides;

    // Utilities
    static vector<int> inferShape(const vector<double>& data);
    static vector<int> inferShape(const vector<vector<double>>& data);
    static vector<int> inferShape(const vector<vector<vector<double>>>& data);

    static vector<int> computeStrides(const vector<int>& shape);
    static int computeNumberOfElements(const vector<int>& shape);
    void validateIndices(const vector<int>& indices) const;
    static vector<int> unflattenIndex(int flatIndex, const vector<int>& shape);

    // Broadcasting
    static vector<int> broadcastShape(const vector<int>& shape1, const vector<int>& shape2);

    // Internal for broadcast_to
    [[nodiscard]] double getBroadcasted(const vector<int>& indices, const vector<int>& expanded_shape) const;
    bool operator<(const Tensor &tensor) const{
        for (int i = 0; i < this->data.size() && i < tensor.data.size(); i++) {
            if (this->data[i] < tensor.data[i]) {
                return true;
            }
            if (this->data[i] > tensor.data[i]) {
                return false;
            }
        }
        return this->data.size() < tensor.data.size();
    }
};

#endif //MATH_TENSOR_H
