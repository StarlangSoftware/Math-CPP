//
// Created by Yiğit Demirşan on 21.02.2025.
//

#ifndef MATH_TENSOR_H
#define MATH_TENSOR_H

#include <initializer_list>
#include <string>

using namespace std;

/**
 * A multidimensional tensor which stores its data, shape and strides in primitive arrays.
 */
class Tensor {
public:
    Tensor();
    /**
     * Takes the ownership of the given data array, therefore the data array must not be deleted or reused
     * by the caller. The shape array is copied.
     */
    Tensor(double* data, const int* shape, int shapeSize);
    Tensor(double* data, initializer_list<int> shape);
    Tensor(initializer_list<double> data, initializer_list<int> shape);
    explicit Tensor(initializer_list<double> data);
    explicit Tensor(initializer_list<initializer_list<double>> data);
    explicit Tensor(initializer_list<initializer_list<initializer_list<double>>> data);
    Tensor(const Tensor& tensor);
    Tensor(Tensor&& tensor) noexcept;
    Tensor& operator=(const Tensor& tensor);
    Tensor& operator=(Tensor&& tensor) noexcept;
    ~Tensor();

    // Element access
    [[nodiscard]] double getValue(const int* indices) const;
    [[nodiscard]] double getValue(initializer_list<int> indices) const;
    void setValue(const int* indices, double value);
    void setValue(initializer_list<int> indices, double value);

    // Data, shape and size
    [[nodiscard]] double* getData() const;
    [[nodiscard]] int* getShape() const;
    [[nodiscard]] int getShapeSize() const;
    [[nodiscard]] int getSize() const;

    // Reshape and transpose
    [[nodiscard]] Tensor reshape(const int* newShape, int newShapeSize) const;
    [[nodiscard]] Tensor reshape(initializer_list<int> newShape) const;
    [[nodiscard]] Tensor transpose(const int* axes = nullptr) const;
    [[nodiscard]] Tensor transpose(initializer_list<int> axes) const;

    [[nodiscard]] Tensor concat(const Tensor& tensor, int dimension) const;
    [[nodiscard]] Tensor get(const int* dimensions, int dimensionSize) const;
    [[nodiscard]] Tensor get(initializer_list<int> dimensions) const;

    // Broadcasting and elementwise ops
    [[nodiscard]] Tensor broadcastTo(const int* targetShape, int targetShapeSize) const;
    [[nodiscard]] Tensor broadcastTo(initializer_list<int> targetShape) const;
    [[nodiscard]] Tensor add(const Tensor& other) const;
    [[nodiscard]] Tensor subtract(const Tensor& other) const;
    [[nodiscard]] Tensor hadamardProduct(const Tensor& other) const;

    // Matrix multiplication
    [[nodiscard]] Tensor multiply(const Tensor& other) const;

    // Slicing
    [[nodiscard]] Tensor partial(const int* startIndices, const int* endIndices) const;
    [[nodiscard]] Tensor partial(initializer_list<int> startIndices, initializer_list<int> endIndices) const;

    // Representation
    [[nodiscard]] string to_string() const;
    bool operator<(const Tensor& tensor) const;
private:
    double* data;
    int* shape;
    int* strides;
    int shapeSize;
    int size;

    void initialize(double* data, const int* shape, int shapeSize);

    // Utilities
    static int* copyOf(const int* array, int size);
    static int* copyOf(initializer_list<int> array);
    static int computeNumberOfElements(const int* shape, int shapeSize);
    static int* computeStrides(const int* shape, int shapeSize);
    static void unflattenIndex(int flatIndex, const int* strides, int stridesSize, int* indices);
    void validateIndices(const int* indices, int indicesSize) const;
    [[nodiscard]] int flatIndex(const int* indices) const;

    static string shapeToString(const int* shape, int shapeSize);
    static string doubleToString(double value);
    static string formatTensor(const double* dataArray, const int* shape, int shapeSize, int offset);

    // Broadcasting
    static bool equalShapes(const int* shape1, int shapeSize1, const int* shape2, int shapeSize2);
    [[nodiscard]] Tensor elementWise(const Tensor& other, double (*operation)(double, double)) const;
    static int* broadcastShape(const int* shape1, int shapeSize1, const int* shape2, int shapeSize2);
    [[nodiscard]] double getBroadcasted(const int* indices, int expandedShapeSize) const;
};

#endif //MATH_TENSOR_H
