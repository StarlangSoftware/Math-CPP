//
// Created by Yiğit Demirşan on 21.02.2025.
//

#include "Tensor.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace std;

/**
 * Initializes an empty tensor.
 */
Tensor::Tensor() {
    this->data = nullptr;
    this->shape = nullptr;
    this->strides = nullptr;
    this->shapeSize = 0;
    this->size = 0;
}

/**
 * Sets the data, shape and strides of the tensor. The ownership of the data array is taken, the shape array
 * is copied.
 *
 * @param data Flattened array representing the tensor data.
 * @param shape The shape of the tensor.
 * @param shapeSize Number of dimensions of the tensor.
 */
void Tensor::initialize(double *data, const int *shape, const int shapeSize) {
    this->data = data;
    this->shapeSize = shapeSize;
    this->shape = copyOf(shape, shapeSize);
    this->strides = computeStrides(this->shape, shapeSize);
    this->size = computeNumberOfElements(this->shape, shapeSize);
}

/**
 * Initializes the tensor with the given data and shape. The tensor takes the ownership of the data array,
 * therefore the caller must neither delete nor reuse it.
 *
 * @param data Flattened array representing the tensor data.
 * @param shape The shape of the tensor.
 * @param shapeSize Number of dimensions of the tensor.
 */
Tensor::Tensor(double* data, const int* shape, const int shapeSize) {
    initialize(data, shape, shapeSize);
}

/**
 * Initializes the tensor with the given data and shape. The tensor takes the ownership of the data array.
 *
 * @param data Flattened array representing the tensor data.
 * @param shape The shape of the tensor.
 */
Tensor::Tensor(double *data, const initializer_list<int> shape) {
    initialize(data, shape.begin(), static_cast<int>(shape.size()));
}

/**
 * Initializes the tensor by copying the given data and shape.
 *
 * @param data Flattened data of the tensor.
 * @param shape The shape of the tensor.
 */
Tensor::Tensor(const initializer_list<double> data, const initializer_list<int> shape) {
    if (computeNumberOfElements(shape.begin(), static_cast<int>(shape.size())) != static_cast<int>(data.size())) {
        throw invalid_argument("Shape does not match the number of elements in data.");
    }
    auto* newData = new double[data.size()];
    copy(data.begin(), data.end(), newData);
    initialize(newData, shape.begin(), static_cast<int>(shape.size()));
}

/**
 * Initializes a one dimensional tensor by copying the given data.
 *
 * @param data Data of the tensor.
 */
Tensor::Tensor(const initializer_list<double> data) {
    auto* newData = new double[data.size()];
    copy(data.begin(), data.end(), newData);
    const int newShape[1] = {static_cast<int>(data.size())};
    initialize(newData, newShape, 1);
}

/**
 * Initializes a two dimensional tensor by flattening the given nested lists.
 *
 * @param data Nested lists representing the tensor data.
 */
Tensor::Tensor(const initializer_list<initializer_list<double>> data) {
    const int rowCount = static_cast<int>(data.size());
    const int columnCount = rowCount == 0 ? 0 : static_cast<int>(data.begin()->size());
    auto* newData = new double[rowCount * columnCount];
    int index = 0;
    for (const auto& row : data) {
        if (static_cast<int>(row.size()) != columnCount) {
            delete[] newData;
            throw invalid_argument("Shape does not match the number of elements in data.");
        }
        for (const double value : row) {
            newData[index++] = value;
        }
    }
    const int newShape[2] = {rowCount, columnCount};
    initialize(newData, newShape, 2);
}

/**
 * Initializes a three dimensional tensor by flattening the given nested lists.
 *
 * @param data Nested lists representing the tensor data.
 */
Tensor::Tensor(const initializer_list<initializer_list<initializer_list<double>>> data) {
    const int firstCount = static_cast<int>(data.size());
    const int secondCount = firstCount == 0 ? 0 : static_cast<int>(data.begin()->size());
    const int thirdCount = secondCount == 0 ? 0 : static_cast<int>(data.begin()->begin()->size());
    auto* newData = new double[firstCount * secondCount * thirdCount];
    int index = 0;
    for (const auto& matrix : data) {
        if (static_cast<int>(matrix.size()) != secondCount) {
            delete[] newData;
            throw invalid_argument("Shape does not match the number of elements in data.");
        }
        for (const auto& row : matrix) {
            if (static_cast<int>(row.size()) != thirdCount) {
                delete[] newData;
                throw invalid_argument("Shape does not match the number of elements in data.");
            }
            for (const double value : row) {
                newData[index++] = value;
            }
        }
    }
    const int newShape[3] = {firstCount, secondCount, thirdCount};
    initialize(newData, newShape, 3);
}

Tensor::Tensor(const Tensor &tensor) {
    this->shapeSize = tensor.shapeSize;
    this->size = tensor.size;
    if (tensor.data == nullptr) {
        this->data = nullptr;
        this->shape = nullptr;
        this->strides = nullptr;
    } else {
        this->data = new double[tensor.size];
        copy(tensor.data, tensor.data + tensor.size, this->data);
        this->shape = copyOf(tensor.shape, tensor.shapeSize);
        this->strides = copyOf(tensor.strides, tensor.shapeSize);
    }
}

Tensor::Tensor(Tensor &&tensor) noexcept {
    this->data = tensor.data;
    this->shape = tensor.shape;
    this->strides = tensor.strides;
    this->shapeSize = tensor.shapeSize;
    this->size = tensor.size;
    tensor.data = nullptr;
    tensor.shape = nullptr;
    tensor.strides = nullptr;
    tensor.shapeSize = 0;
    tensor.size = 0;
}

Tensor & Tensor::operator=(const Tensor &tensor) {
    if (this == &tensor) {
        return *this;
    }
    delete[] this->data;
    delete[] this->shape;
    delete[] this->strides;
    this->shapeSize = tensor.shapeSize;
    this->size = tensor.size;
    if (tensor.data == nullptr) {
        this->data = nullptr;
        this->shape = nullptr;
        this->strides = nullptr;
    } else {
        this->data = new double[tensor.size];
        copy(tensor.data, tensor.data + tensor.size, this->data);
        this->shape = copyOf(tensor.shape, tensor.shapeSize);
        this->strides = copyOf(tensor.strides, tensor.shapeSize);
    }
    return *this;
}

Tensor & Tensor::operator=(Tensor &&tensor) noexcept {
    if (this == &tensor) {
        return *this;
    }
    delete[] this->data;
    delete[] this->shape;
    delete[] this->strides;
    this->data = tensor.data;
    this->shape = tensor.shape;
    this->strides = tensor.strides;
    this->shapeSize = tensor.shapeSize;
    this->size = tensor.size;
    tensor.data = nullptr;
    tensor.shape = nullptr;
    tensor.strides = nullptr;
    tensor.shapeSize = 0;
    tensor.size = 0;
    return *this;
}

Tensor::~Tensor() {
    delete[] this->data;
    delete[] this->shape;
    delete[] this->strides;
}

// --- Data and shape ---

double * Tensor::getData() const {
    return this->data;
}

int * Tensor::getShape() const {
    return this->shape;
}

int Tensor::getShapeSize() const {
    return this->shapeSize;
}

int Tensor::getSize() const {
    return this->size;
}

/**
 * Copies the given array into a newly allocated array.
 *
 * @param array Array to be copied.
 * @param size Size of the array.
 * @return The newly allocated copy of the array.
 */
int * Tensor::copyOf(const int *array, const int size) {
    if (array == nullptr) {
        return nullptr;
    }
    auto* newArray = new int[size];
    copy(array, array + size, newArray);
    return newArray;
}

int * Tensor::copyOf(const initializer_list<int> array) {
    auto* newArray = new int[array.size()];
    copy(array.begin(), array.end(), newArray);
    return newArray;
}

/**
 * Computes the strides for each dimension based on the shape.
 *
 * @param shape Array representing the tensor shape.
 * @param shapeSize Number of dimensions.
 * @return Array representing the strides.
 */
int* Tensor::computeStrides(const int* shape, const int shapeSize) {
    auto* strides = new int[shapeSize];
    int prod = 1;
    for (int i = shapeSize - 1; i >= 0; --i) {
        strides[i] = prod;
        prod *= shape[i];
    }
    return strides;
}

/**
 * Computes the total number of elements in the tensor based on its shape.
 *
 * @param shape Array representing the tensor shape.
 * @param shapeSize Number of dimensions.
 * @return Total number of elements.
 */
int Tensor::computeNumberOfElements(const int* shape, const int shapeSize) {
    int product = 1;
    for (int i = 0; i < shapeSize; i++) {
        product *= shape[i];
    }
    return product;
}

/**
 * Returns the string representation of the given shape such as [2, 3].
 *
 * @param shape Array representing the tensor shape.
 * @param shapeSize Number of dimensions.
 * @return String representation of the shape.
 */
string Tensor::shapeToString(const int* shape, const int shapeSize) {
    ostringstream oss;
    oss << "[";
    for (int i = 0; i < shapeSize; ++i) {
        oss << shape[i];
        if (i + 1 < shapeSize) oss << ", ";
    }
    oss << "]";
    return oss.str();
}

/**
 * Validates that indices are within the valid range for each dimension.
 *
 * @param indices Array of indices specifying the position.
 * @param indicesSize Number of indices.
 */
void Tensor::validateIndices(const int* indices, const int indicesSize) const {
    if (indicesSize != shapeSize)
        throw out_of_range("Expected " + std::to_string(shapeSize) + " indices but got " + std::to_string(indicesSize) + ".");
    for (int i = 0; i < indicesSize; ++i) {
        if (indices[i] < 0 || indices[i] >= shape[i])
            throw out_of_range("Index " + shapeToString(indices, indicesSize) + " is out of bounds for shape " + shapeToString(shape, shapeSize) + ".");
    }
}

/**
 * Retrieves the value at the given indices.
 *
 * @param indices Array of indices specifying the position. It must contain one index for each dimension.
 * @return Value at the specified position.
 */
int Tensor::flatIndex(const int *indices) const {
    int flat = 0;
    for (int i = 0; i < shapeSize; ++i) {
        flat += indices[i] * strides[i];
    }
    return flat;
}

double Tensor::getValue(const int* indices) const {
    validateIndices(indices, shapeSize);
    int flat = 0;
    for (int i = 0; i < shapeSize; ++i) {
        flat += indices[i] * strides[i];
    }
    return data[flat];
}

double Tensor::getValue(const initializer_list<int> indices) const {
    validateIndices(indices.begin(), static_cast<int>(indices.size()));
    int flat = 0;
    int i = 0;
    for (const int index : indices) {
        flat += index * strides[i++];
    }
    return data[flat];
}

/**
 * Sets the value at the given indices.
 *
 * @param indices Array of indices specifying the position. It must contain one index for each dimension.
 * @param value   Value to set at the specified position.
 */
void Tensor::setValue(const int* indices, const double value) {
    validateIndices(indices, shapeSize);
    int flat = 0;
    for (int i = 0; i < shapeSize; ++i) {
        flat += indices[i] * strides[i];
    }
    data[flat] = value;
}

void Tensor::setValue(const initializer_list<int> indices, const double value) {
    validateIndices(indices.begin(), static_cast<int>(indices.size()));
    int flat = 0;
    int i = 0;
    for (const int index : indices) {
        flat += index * strides[i++];
    }
    data[flat] = value;
}

/**
 * Converts a flat index to multidimensional indices based on the given strides.
 *
 * @param flatIndex The flat index to convert.
 * @param strides Array representing the strides.
 * @param stridesSize Number of dimensions.
 * @param indices Array in which the multi-dimensional indices are written.
 */
void Tensor::unflattenIndex(int flatIndex, const int* strides, const int stridesSize, int* indices) {
    for (int i = 0; i < stridesSize; ++i) {
        indices[i] = flatIndex / strides[i];
        flatIndex %= strides[i];
    }
}

/**
 * Reshapes the tensor to the specified new shape.
 *
 * @param newShape Array representing the new shape.
 * @param newShapeSize Number of dimensions of the new shape.
 * @return New tensor with the specified shape.
 */
Tensor Tensor::reshape(const int* newShape, const int newShapeSize) const {
    if (computeNumberOfElements(newShape, newShapeSize) != size)
        throw invalid_argument("Total number of elements must remain the same.");
    auto* newData = new double[size];
    copy(data, data + size, newData);
    return {newData, newShape, newShapeSize};
}

Tensor Tensor::reshape(const initializer_list<int> newShape) const {
    return reshape(newShape.begin(), static_cast<int>(newShape.size()));
}

/**
 * Transposes the tensor according to the specified axes.
 *
 * @param axes Array representing the order of axes. If nullptr, reverses the axes.
 * @return New tensor with transposed axes.
 */
Tensor Tensor::transpose(const int* axes) const {
    auto* axesOrder = new int[shapeSize];
    if (axes == nullptr) {
        for (int i = 0; i < shapeSize; ++i) {
            axesOrder[i] = shapeSize - 1 - i;
        }
    } else {
        copy(axes, axes + shapeSize, axesOrder);
        auto* sortedAxes = copyOf(axesOrder, shapeSize);
        sort(sortedAxes, sortedAxes + shapeSize);
        for (int i = 0; i < shapeSize; ++i) {
            if (sortedAxes[i] != i) {
                delete[] sortedAxes;
                delete[] axesOrder;
                throw invalid_argument("Invalid transpose axes.");
            }
        }
        delete[] sortedAxes;
    }
    auto* newShape = new int[shapeSize];
    for (int i = 0; i < shapeSize; ++i) {
        newShape[i] = shape[axesOrder[i]];
    }
    int* newStrides = computeStrides(newShape, shapeSize);
    auto* newData = new double[size];
    auto* newIndices = new int[shapeSize];
    auto* originalIndices = new int[shapeSize];
    for (int i = 0; i < size; ++i) {
        unflattenIndex(i, newStrides, shapeSize, newIndices);
        for (int j = 0; j < shapeSize; ++j) {
            originalIndices[axesOrder[j]] = newIndices[j];
        }
        newData[i] = getValue(originalIndices);
    }
    delete[] newIndices;
    delete[] originalIndices;
    delete[] newStrides;
    delete[] axesOrder;
    Tensor result = Tensor(newData, newShape, shapeSize);
    delete[] newShape;
    return result;
}

Tensor Tensor::transpose(const initializer_list<int> axes) const {
    if (static_cast<int>(axes.size()) != shapeSize) {
        throw invalid_argument("Invalid transpose axes.");
    }
    return transpose(axes.begin());
}

/**
 * Concatenates two tensors into a one.
 *
 * @param tensor 2nd tensor for concatenation.
 * @param dimension to concatenate.
 * @return Concatenated {@link Tensor}.
 */
Tensor Tensor::concat(const Tensor &tensor, const int dimension) const {
    if (dimension >= shapeSize) {
        throw invalid_argument("Dimension out of bounds.");
    }
    if (shapeSize != tensor.shapeSize) {
        throw invalid_argument("Dimensions length do not match.");
    }
    int startIndex = 1;
    int endIndex1 = 1;
    int endIndex2 = 1;
    for (int i = 0; i < shapeSize; i++) {
        if (i != dimension && shape[i] != tensor.shape[i]) {
            throw invalid_argument("Dimensions do not match.");
        }
        if (i >= dimension) {
            endIndex1 *= shape[i];
            endIndex2 *= tensor.shape[i];
        } else {
            startIndex *= shape[i];
        }
    }
    auto* newShape = new int[shapeSize];
    for (int i = 0; i < shapeSize; i++) {
        if (i == dimension) {
            newShape[i] = shape[i] + tensor.shape[i];
        } else {
            newShape[i] = shape[i];
        }
    }
    auto* newData = new double[startIndex * (endIndex1 + endIndex2)];
    int k = 0;
    for (int i = 0; i < startIndex; i++) {
        for (int j = 0; j < endIndex1; j++) {
            newData[k++] = data[i * endIndex1 + j];
        }
        for (int j = 0; j < endIndex2; j++) {
            newData[k++] = tensor.data[i * endIndex2 + j];
        }
    }
    Tensor result = Tensor(newData, newShape, shapeSize);
    delete[] newShape;
    return result;
}

/**
 * Returns the sub-{@link Tensor} taking the given dimensions.
 *
 * @param dimensions Indices of the dimensions to be taken.
 * @param dimensionSize Number of the given dimensions.
 * @return a sub-{@link Tensor}.
 */
Tensor Tensor::get(const int* dimensions, const int dimensionSize) const {
    if (dimensionSize >= shapeSize) {
        throw invalid_argument("Dimensions exceeds or same as the tensor's dimension.");
    }
    for (int i = 0; i < dimensionSize; i++) {
        if (dimensions[i] >= shape[i]) {
            throw out_of_range("There is a dimension length exceed the tensor's dimension length.");
        }
    }
    int i = 0, start = 0, end = size;
    do {
        const int parts = (end - start) / shape[i];
        start += parts * dimensions[i];
        end = start + parts;
        i++;
    } while (i < dimensionSize);
    auto* newData = new double[end - start];
    copy(data + start, data + end, newData);
    return {newData, shape + dimensionSize, shapeSize - dimensionSize};
}

Tensor Tensor::get(const initializer_list<int> dimensions) const {
    return get(dimensions.begin(), static_cast<int>(dimensions.size()));
}

/**
 * Determines the broadcasted shape of two tensors.
 *
 * @param shape1 Array representing the first tensor shape.
 * @param shapeSize1 Number of dimensions of the first tensor.
 * @param shape2 Array representing the second tensor shape.
 * @param shapeSize2 Number of dimensions of the second tensor.
 * @return Array representing the broadcasted shape. Its size is the maximum of the given dimension counts.
 */
int* Tensor::broadcastShape(const int* shape1, const int shapeSize1, const int* shape2, const int shapeSize2) {
    const int ndim = max(shapeSize1, shapeSize2);
    const int offset1 = ndim - shapeSize1;
    const int offset2 = ndim - shapeSize2;
    auto* out = new int[ndim];
    for (int i = 0; i < ndim; ++i) {
        const int dimension1 = i < offset1 ? 1 : shape1[i - offset1];
        const int dimension2 = i < offset2 ? 1 : shape2[i - offset2];
        if (dimension1 == dimension2 || dimension2 == 1) {
            out[i] = dimension1;
        } else if (dimension1 == 1) {
            out[i] = dimension2;
        } else {
            delete[] out;
            throw invalid_argument("Shapes " + shapeToString(shape1, shapeSize1) + " and " + shapeToString(shape2, shapeSize2) + " are not broadcastable");
        }
    }
    return out;
}

double Tensor::getBroadcasted(const int* indices, const int expandedShapeSize) const {
    const int dimensionOffset = expandedShapeSize - shapeSize;
    int flat = 0;
    for (int i = 0; i < shapeSize; ++i) {
        if (shape[i] != 1) {
            flat += indices[dimensionOffset + i] * strides[i];
        }
    }
    return data[flat];
}

/**
 * Broadcasts the tensor to the specified target shape.
 *
 * @param targetShape Array representing the target shape.
 * @param targetShapeSize Number of dimensions of the target shape.
 * @return New tensor with the target shape.
 */
Tensor Tensor::broadcastTo(const int* targetShape, const int targetShapeSize) const {
    if (targetShapeSize < shapeSize) {
        throw invalid_argument("Cannot broadcast shape " + shapeToString(shape, shapeSize) + " to " + shapeToString(targetShape, targetShapeSize));
    }
    const int dimensionOffset = targetShapeSize - shapeSize;
    for (int i = 0; i < targetShapeSize; ++i) {
        const int expandedDimension = i < dimensionOffset ? 1 : shape[i - dimensionOffset];
        if (expandedDimension != targetShape[i] && expandedDimension != 1) {
            throw invalid_argument("Cannot broadcast shape " + shapeToString(shape, shapeSize) + " to " + shapeToString(targetShape, targetShapeSize));
        }
    }
    const int total = computeNumberOfElements(targetShape, targetShapeSize);
    int* targetStrides = computeStrides(targetShape, targetShapeSize);
    auto* newData = new double[total];
    auto* indices = new int[targetShapeSize];
    for (int idx = 0; idx < total; ++idx) {
        unflattenIndex(idx, targetStrides, targetShapeSize, indices);
        newData[idx] = getBroadcasted(indices, targetShapeSize);
    }
    delete[] indices;
    delete[] targetStrides;
    return {newData, targetShape, targetShapeSize};
}

Tensor Tensor::broadcastTo(const initializer_list<int> targetShape) const {
    return broadcastTo(targetShape.begin(), static_cast<int>(targetShape.size()));
}

/**
 * Checks whether the two shapes are identical.
 *
 * @param shape1 Array representing the first shape.
 * @param shapeSize1 Number of dimensions of the first shape.
 * @param shape2 Array representing the second shape.
 * @param shapeSize2 Number of dimensions of the second shape.
 * @return True if the shapes are identical, false otherwise.
 */
bool Tensor::equalShapes(const int *shape1, const int shapeSize1, const int *shape2, const int shapeSize2) {
    if (shapeSize1 != shapeSize2) {
        return false;
    }
    for (int i = 0; i < shapeSize1; i++) {
        if (shape1[i] != shape2[i]) {
            return false;
        }
    }
    return true;
}

/**
 * Applies the given element-wise operation to the two tensors, broadcasting them only if their shapes differ.
 *
 * @param other The other tensor.
 * @param operation Element-wise operation to be applied.
 * @return New tensor with the result of the operation.
 */
Tensor Tensor::elementWise(const Tensor &other, double (*operation)(double, double)) const {
    if (equalShapes(shape, shapeSize, other.shape, other.shapeSize)) {
        auto* result = new double[size];
        for (int i = 0; i < size; ++i) {
            result[i] = operation(data[i], other.data[i]);
        }
        return {result, shape, shapeSize};
    }
    const int outShapeSize = max(shapeSize, other.shapeSize);
    int* outShape = broadcastShape(shape, shapeSize, other.shape, other.shapeSize);
    const Tensor first = broadcastTo(outShape, outShapeSize);
    const Tensor second = other.broadcastTo(outShape, outShapeSize);
    auto* result = new double[first.size];
    for (int i = 0; i < first.size; ++i) {
        result[i] = operation(first.data[i], second.data[i]);
    }
    Tensor resultTensor = Tensor(result, outShape, outShapeSize);
    delete[] outShape;
    return resultTensor;
}

/**
 * Adds two tensors element-wise with broadcasting.
 *
 * @param other The other tensor to add.
 * @return New tensor with the result of the addition.
 */
Tensor Tensor::add(const Tensor& other) const {
    return elementWise(other, [](const double first, const double second) { return first + second; });
}

/**
 * Subtracts one tensor from another element-wise with broadcasting.
 *
 * @param other The other tensor to subtract.
 * @return New tensor with the result of the subtraction.
 */
Tensor Tensor::subtract(const Tensor& other) const {
    return elementWise(other, [](const double first, const double second) { return first - second; });
}

/**
 * Multiplies two tensors element-wise with broadcasting.
 *
 * @param other The other tensor to multiply.
 * @return New tensor with the result of the multiplication.
 */
Tensor Tensor::hadamardProduct(const Tensor& other) const {
    return elementWise(other, [](const double first, const double second) { return first * second; });
}

/**
 * Performs matrix multiplication (batched if necessary).
 * For tensors of shape (..., M, K) and (..., K, N), returns (..., M, N).
 *
 * @param other Tensor with shape compatible for matrix multiplication.
 * @return Tensor resulting from matrix multiplication.
 */
Tensor Tensor::multiply(const Tensor& other) const {
    if (shapeSize < 2 || other.shapeSize < 2) {
        throw invalid_argument("Shapes " + shapeToString(shape, shapeSize) + " and " + shapeToString(other.shape, other.shapeSize) + " are not aligned for multiplication.");
    }
    if (shape[shapeSize - 1] != other.shape[other.shapeSize - 2]) {
        throw invalid_argument("Shapes " + shapeToString(shape, shapeSize) + " and " + shapeToString(other.shape, other.shapeSize) + " are not aligned for multiplication.");
    }
    const int m = shape[shapeSize - 2];
    const int k = shape[shapeSize - 1];
    const int n = other.shape[other.shapeSize - 1];
    const int batchSize1 = shapeSize - 2;
    const int batchSize2 = other.shapeSize - 2;
    const bool sameBatch = equalShapes(shape, batchSize1, other.shape, batchSize2);
    int batchSize;
    int* batchShape;
    Tensor leftBroadcasted;
    Tensor rightBroadcasted;
    const Tensor* left;
    const Tensor* right;
    if (sameBatch) {
        batchSize = batchSize1;
        batchShape = copyOf(shape, batchSize);
        left = this;
        right = &other;
    } else {
        batchSize = max(batchSize1, batchSize2);
        batchShape = broadcastShape(shape, batchSize1, other.shape, batchSize2);
        auto* leftShape = new int[batchSize + 2];
        auto* rightShape = new int[batchSize + 2];
        copy(batchShape, batchShape + batchSize, leftShape);
        copy(batchShape, batchShape + batchSize, rightShape);
        leftShape[batchSize] = m;
        leftShape[batchSize + 1] = k;
        rightShape[batchSize] = k;
        rightShape[batchSize + 1] = n;
        leftBroadcasted = broadcastTo(leftShape, batchSize + 2);
        rightBroadcasted = other.broadcastTo(rightShape, batchSize + 2);
        left = &leftBroadcasted;
        right = &rightBroadcasted;
        delete[] leftShape;
        delete[] rightShape;
    }
    const int resultShapeSize = batchSize + 2;
    auto* resultShape = new int[resultShapeSize];
    copy(batchShape, batchShape + batchSize, resultShape);
    resultShape[batchSize] = m;
    resultShape[batchSize + 1] = n;
    delete[] batchShape;
    const int resultTotal = computeNumberOfElements(resultShape, resultShapeSize);
    auto* resultData = new double[resultTotal];
    const double* leftData = left->data;
    const double* rightData = right->data;
    const int batchCount = m * n == 0 ? 0 : resultTotal / (m * n);
    for (int batch = 0; batch < batchCount; batch++) {
        const int leftOffset = batch * m * k;
        const int rightOffset = batch * k * n;
        const int resultOffset = batch * m * n;
        for (int row = 0; row < m; row++) {
            for (int col = 0; col < n; col++) {
                double sum = 0.0;
                for (int index = 0; index < k; index++) {
                    sum += leftData[leftOffset + row * k + index] * rightData[rightOffset + index * n + col];
                }
                resultData[resultOffset + row * n + col] = sum;
            }
        }
    }
    Tensor result = Tensor(resultData, resultShape, resultShapeSize);
    delete[] resultShape;
    return result;
}

/**
 * Extracts a sub-tensor from the given start indices to the end indices.
 *
 * @param startIndices Array specifying the start indices for each dimension.
 * @param endIndices   Array specifying the end indices (exclusive) for each dimension.
 * @return A new Tensor containing the extracted sub-tensor.
 */
Tensor Tensor::partial(const int* startIndices, const int* endIndices) const {
    auto* newShape = new int[shapeSize];
    for (int i = 0; i < shapeSize; ++i) {
        newShape[i] = endIndices[i] - startIndices[i];
        if (newShape[i] < 0) {
            delete[] newShape;
            throw invalid_argument("End index must be greater than or equal to start index for dimension " + std::to_string(i));
        }
    }
    for (int i = 0; i < shapeSize; ++i) {
        if (startIndices[i] < 0 || endIndices[i] > shape[i]) {
            const int invalidIndices[1] = {endIndices[i] > shape[i] ? endIndices[i] - 1 : startIndices[i]};
            delete[] newShape;
            throw out_of_range("Index " + shapeToString(invalidIndices, 1) + " is out of bounds for shape " + shapeToString(shape, shapeSize) + ".");
        }
    }
    const int total = computeNumberOfElements(newShape, shapeSize);
    int* newStrides = computeStrides(newShape, shapeSize);
    auto* newData = new double[total];
    auto* subIndices = new int[shapeSize];
    auto* originalIndices = new int[shapeSize];
    for (int idx = 0; idx < total; ++idx) {
        unflattenIndex(idx, newStrides, shapeSize, subIndices);
        for (int j = 0; j < shapeSize; ++j) {
            originalIndices[j] = startIndices[j] + subIndices[j];
        }
        newData[idx] = data[flatIndex(originalIndices)];
    }
    delete[] subIndices;
    delete[] originalIndices;
    delete[] newStrides;
    Tensor result = Tensor(newData, newShape, shapeSize);
    delete[] newShape;
    return result;
}

Tensor Tensor::partial(const initializer_list<int> startIndices, const initializer_list<int> endIndices) const {
    if (static_cast<int>(startIndices.size()) != shapeSize || static_cast<int>(endIndices.size()) != shapeSize) {
        throw invalid_argument("startIndices and endIndices must match the number of dimensions.");
    }
    return partial(startIndices.begin(), endIndices.begin());
}

/**
 * Returns the string representation of a double value in the same format with Java.
 *
 * @param value Value to be converted.
 * @return String representation of the value.
 */
string Tensor::doubleToString(const double value) {
    if (isnan(value)) {
        return "NaN";
    }
    if (isinf(value)) {
        return value > 0 ? "Infinity" : "-Infinity";
    }
    string result;
    for (int precision = 1; precision <= 17; precision++) {
        ostringstream oss;
        oss << setprecision(precision) << value;
        result = oss.str();
        if (stod(result) == value) {
            break;
        }
    }
    if (result.find('.') == string::npos && result.find('e') == string::npos && result.find('E') == string::npos) {
        result += ".0";
    }
    return result;
}

/**
 * Formats the data of the tensor as nested lists.
 *
 * @param dataArray Flattened data of the tensor.
 * @param shape Shape of the (sub-)tensor to be formatted.
 * @param shapeSize Number of dimensions of the (sub-)tensor.
 * @param offset Start index of the (sub-)tensor in the flattened data.
 * @return String representation of the data.
 */
string Tensor::formatTensor(const double* dataArray, const int* shape, const int shapeSize, const int offset) {
    if (shapeSize == 0) {
        return "[]";
    }
    ostringstream oss;
    oss << "[";
    if (shapeSize == 1) {
        for (int i = 0; i < shape[0]; ++i) {
            oss << doubleToString(dataArray[offset + i]);
            if (i + 1 < shape[0]) oss << ", ";
        }
    } else {
        const int stride = computeNumberOfElements(shape + 1, shapeSize - 1);
        for (int i = 0; i < shape[0]; ++i) {
            oss << formatTensor(dataArray, shape + 1, shapeSize - 1, offset + (i * stride));
            if (i + 1 < shape[0]) oss << ", ";
        }
    }
    oss << "]";
    return oss.str();
}

/**
 * Returns a string representation of the tensor.
 *
 * @return String representing the tensor.
 */
string Tensor::to_string() const {
    return "Tensor(shape=" + shapeToString(shape, shapeSize) + ", data=" + formatTensor(data, shape, shapeSize, 0) + ")";
}

bool Tensor::operator<(const Tensor &tensor) const {
    for (int i = 0; i < size && i < tensor.size; i++) {
        if (data[i] < tensor.data[i]) {
            return true;
        }
        if (data[i] > tensor.data[i]) {
            return false;
        }
    }
    return size < tensor.size;
}
