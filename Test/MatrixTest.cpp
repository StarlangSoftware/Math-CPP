//
// Created by Olcay Taner YILDIZ on 23.12.2020.
//

#include "catch.hpp"
#include "../src/Matrix.h"

TEST_CASE("MatrixTest-testColumnWiseNormalize") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++){
        for (int j = 0; j < 1000; j++){
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    small.columnWiseNormalize();
    REQUIRE(3 == small.sumOfElements());
    large.columnWiseNormalize();
    REQUIRE_THAT(1000, Catch::Matchers::WithinAbs(large.sumOfElements(), 0.001));
    identity.columnWiseNormalize();
    REQUIRE(100 == identity.sumOfElements());
}

TEST_CASE("MatrixTest-testMultiplyWithConstant") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    small.multiplyWithConstant(4);
    REQUIRE(36 == small.sumOfElements());
    large.multiplyWithConstant(1.001);
    REQUIRE_THAT(1001000, Catch::Matchers::WithinAbs(large.sumOfElements(), 0.001));
    random.multiplyWithConstant(3.6);
    REQUIRE_THAT(originalSum * 3.6, Catch::Matchers::WithinAbs(random.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testDivideByConstant") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    small.divideByConstant(4);
    REQUIRE(2.25 == small.sumOfElements());
    large.divideByConstant(10);
    REQUIRE_THAT(100000, Catch::Matchers::WithinAbs(large.sumOfElements(), 0.001));
    random.divideByConstant(3.6);
    REQUIRE_THAT(originalSum / 3.6, Catch::Matchers::WithinAbs(random.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testAdd") {
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    auto identity = Matrix(100);
    random.add(identity);
    REQUIRE_THAT(originalSum + 100, Catch::Matchers::WithinAbs(random.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testAddVector") {
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    Vector V = Vector(1000, 1.0);
    large.add(4, V);
    REQUIRE(1001000 == large.sumOfElements());
}

TEST_CASE("MatrixTest-testSubtract") {
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    auto identity = Matrix(100);
    random.subtract(identity);
    REQUIRE_THAT(originalSum - 100, Catch::Matchers::WithinAbs(random.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testMultiplyWithVectorFromLeft") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    Vector v = Vector(3, 1.0);
    Vector V = Vector(1000, 1.0);
    Vector vr = Vector(100, 1.0);
    Vector result = small.multiplyWithVectorFromLeft(v);
    REQUIRE(9 == result.sumOfElements());
    result = large.multiplyWithVectorFromLeft(V);
    REQUIRE(1000000 == result.sumOfElements());
    result = random.multiplyWithVectorFromLeft(vr);
    REQUIRE_THAT(originalSum, Catch::Matchers::WithinAbs(result.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testMultiplyWithVectorFromRight") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    Vector v = Vector(3, 1.0);
    Vector V = Vector(1000, 1.0);
    Vector vr = Vector(100, 1.0);
    Vector result = small.multiplyWithVectorFromRight(v);
    REQUIRE(9 == result.sumOfElements());
    result = large.multiplyWithVectorFromRight(V);
    REQUIRE(1000000 == result.sumOfElements());
    result = random.multiplyWithVectorFromRight(vr);
    REQUIRE_THAT(originalSum, Catch::Matchers::WithinAbs(result.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testColumnSum") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    REQUIRE(3 == small.columnSum(random() % 3));
    REQUIRE(1000 == large.columnSum(random() % 1000));
    REQUIRE(1 == identity.columnSum(random() % 100));
}

TEST_CASE("MatrixTest-testSumOfRows") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    REQUIRE(9 == small.sumOfRows().sumOfElements());
    REQUIRE(1000000 == large.sumOfRows().sumOfElements());
    REQUIRE(100 == identity.sumOfRows().sumOfElements());
    REQUIRE_THAT(originalSum, Catch::Matchers::WithinAbs(random.sumOfRows().sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testRowSum") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    REQUIRE(3 == small.rowSum(random() % 3));
    REQUIRE(1000 == large.rowSum(random() % 1000));
    REQUIRE(1 == identity.rowSum(random() % 100));
}

TEST_CASE("MatrixTest-testMultiply") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    Matrix result = small.multiply(small);
    REQUIRE(27 == result.sumOfElements());
    Matrix result2 = large.multiply(large);
    REQUIRE(1000000000.0 == result2.sumOfElements());
    Matrix result3 = random.multiply(identity);
    REQUIRE(originalSum == result3.sumOfElements());
    Matrix result4 = identity.multiply(random);
    REQUIRE(originalSum == result4.sumOfElements());
}

TEST_CASE("MatrixTest-testElementProduct") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    Matrix result = small.elementProduct(small);
    REQUIRE(9 == result.sumOfElements());
    Matrix result2 = large.elementProduct(large);
    REQUIRE(1000000 == result2.sumOfElements());
    Matrix result3 = random.elementProduct(identity);
    REQUIRE(result3.trace() == result3.sumOfElements());
}

TEST_CASE("MatrixTest-testSumOfElements") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    REQUIRE(9 == small.sumOfRows().sumOfElements());
    REQUIRE(1000000 == large.sumOfRows().sumOfElements());
    REQUIRE(100 == identity.sumOfRows().sumOfElements());
    REQUIRE_THAT(originalSum, Catch::Matchers::WithinAbs(random.sumOfRows().sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testTrace") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    REQUIRE(3 == small.trace());
    REQUIRE(1000 == large.trace());
    REQUIRE(100 == identity.trace());
}

TEST_CASE("MatrixTest-testTranspose") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    REQUIRE(9 == small.transpose().sumOfElements());
    REQUIRE(1000000 == large.transpose().sumOfElements());
    REQUIRE(100 == identity.transpose().sumOfElements());
    REQUIRE_THAT(originalSum, Catch::Matchers::WithinAbs(random.transpose().sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testIsSymmetric") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    REQUIRE_FALSE(!small.isSymmetric());
    REQUIRE_FALSE(!large.isSymmetric());
    REQUIRE_FALSE(!identity.isSymmetric());
    REQUIRE_FALSE(random.isSymmetric());
}

TEST_CASE("MatrixTest-testDeterminant") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    REQUIRE(0 == small.determinant());
    REQUIRE(0 == large.determinant());
    REQUIRE(1 == identity.determinant());
}

TEST_CASE("MatrixTest-testInverse") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    Matrix large = Matrix(1000, 1000);
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 1000; j++) {
            large.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix random = Matrix(100, 100, 1, 10, default_random_engine(1));
    double originalSum = random.sumOfElements();
    identity.inverse();
    REQUIRE(100 == identity.sumOfElements());
    random.inverse();
    random.inverse();
    REQUIRE_THAT(originalSum, Catch::Matchers::WithinAbs(random.sumOfElements(), 0.001));
}

TEST_CASE("MatrixTest-testCharacteristics") {
    Matrix small = Matrix(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            small.setValue(i, j, 1.0);
        }
    }
    auto identity = Matrix(100);
    Matrix medium = Matrix(100, 100);
    for (int i = 0; i < 100; i++){
        for (int j = 0; j < 100; j++){
            medium.setValue(i, j, 1.0);
        }
    }
    vector<Eigenvector> vectors = small.characteristics();
    REQUIRE(2 == vectors.size());
    vectors = identity.characteristics();
    REQUIRE(100 == vectors.size());
    vectors = medium.characteristics();
    REQUIRE(41 == vectors.size());
}