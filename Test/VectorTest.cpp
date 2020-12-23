//
// Created by Olcay Taner YILDIZ on 23.12.2020.
//

#include "catch.hpp"
#include "../Vector.h"

TEST_CASE("VectorTest-testBiased") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    Vector biased = smallVector1.biased();
    REQUIRE(1 == biased.getValue(0));
    REQUIRE(smallVector1.getSize() + 1 == biased.getSize());
}

TEST_CASE("VectorTest-testElementAdd") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.add(7);
    REQUIRE(7 == smallVector1.getValue(5));
    REQUIRE(6 == smallVector1.getSize());
}

TEST_CASE("VectorTest-testInsert") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.insert(3, 6);
    REQUIRE(6 == smallVector1.getValue(3));
    REQUIRE(6 == smallVector1.getSize());
}

TEST_CASE("VectorTest-testRemove") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.remove(2);
    REQUIRE(5 == smallVector1.getValue(2));
    REQUIRE(4 == smallVector1.getSize());
}

TEST_CASE("VectorTest-testSumOfElementsSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    REQUIRE(20 == smallVector1.sumOfElements());
    REQUIRE(30 == smallVector2.sumOfElements());
}

TEST_CASE("VectorTest-testSumOfElementsLarge") {
    auto* largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++){
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto* largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++){
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    REQUIRE(500500 == largeVector1.sumOfElements());
    REQUIRE(500500 == largeVector2.sumOfElements());
}

TEST_CASE("VectorTest-testMaxIndex") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    REQUIRE(4 == smallVector1.maxIndex());
    REQUIRE(0 == smallVector2.maxIndex());
}

TEST_CASE("VectorTest-testSigmoid") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.sigmoid();
    REQUIRE_THAT(0.8807971, Catch::Matchers::WithinAbs(smallVector1.getValue(0), 0.000001));
    REQUIRE_THAT(0.9975274, Catch::Matchers::WithinAbs(smallVector1.getValue(4), 0.000001));
}

TEST_CASE("VectorTest-testSkipVectorSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector3 = smallVector1.skipVector(2, 0);
    REQUIRE(2 == smallVector3.getValue(0));
    REQUIRE(6 == smallVector3.getValue(2));
    smallVector3 = smallVector1.skipVector(3, 1);
    REQUIRE(3 == smallVector3.getValue(0));
    REQUIRE(6 == smallVector3.getValue(1));
}

TEST_CASE("VectorTest-testSkipVectorLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    Vector largeVector3 = largeVector1.skipVector(2, 0);
    REQUIRE(250000 == largeVector3.sumOfElements());
    largeVector3 = largeVector1.skipVector(5, 3);
    REQUIRE(100300 == largeVector3.sumOfElements());
}

TEST_CASE("VectorTest-testVectorAddSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    smallVector1.add(smallVector2);
    REQUIRE(50 == smallVector1.sumOfElements());
}

TEST_CASE("VectorTest-testVectorAddLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto *largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    largeVector1.add(largeVector2);
    REQUIRE(1001000 == largeVector1.sumOfElements());
}

TEST_CASE("VectorTest-testSubtractSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    smallVector1.subtract(smallVector2);
    REQUIRE(-10 == smallVector1.sumOfElements());
}

TEST_CASE("VectorTest-testSubtractLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto *largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    largeVector1.subtract(largeVector2);
    REQUIRE(0 == largeVector1.sumOfElements());
}

TEST_CASE("VectorTest-testDifferenceSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    Vector smallVector3 = smallVector1.difference(smallVector2);
    REQUIRE(-10 == smallVector3.sumOfElements());
}

TEST_CASE("VectorTest-testDifferenceLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto *largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    Vector largeVector3 = largeVector1.difference(largeVector2);
    REQUIRE(0 == largeVector3.sumOfElements());
}

TEST_CASE("VectorTest-testDotProductWithVectorSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    double dotProduct = smallVector1.dotProduct(smallVector2);
    REQUIRE(110 == dotProduct);
}

TEST_CASE("VectorTest-testDotProductWithVectorLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto *largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    double dotProduct = largeVector1.dotProduct(largeVector2);
    REQUIRE(167167000 == dotProduct);
}

TEST_CASE("VectorTest-testDotProductWithItselfSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    double dotProduct = smallVector1.dotProduct();
    REQUIRE(90 == dotProduct);
}

TEST_CASE("VectorTest-testDotProductWithItselfLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    double dotProduct = largeVector1.dotProduct();
    REQUIRE(333833500 == dotProduct);
}

TEST_CASE("VectorTest-testElementProductSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    Vector smallVector3 = smallVector1.elementProduct(smallVector2);
    REQUIRE(110 == smallVector3.sumOfElements());
}

TEST_CASE("VectorTest-testElementProductLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto *largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    Vector largeVector3 = largeVector1.elementProduct(largeVector2);
    REQUIRE(167167000 == largeVector3.sumOfElements());
}

TEST_CASE("VectorTest-testDivide") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.divide(10.0);
    REQUIRE(2 == smallVector1.sumOfElements());
}

TEST_CASE("VectorTest-testMultiply") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.multiply(10.0);
    REQUIRE(200 == smallVector1.sumOfElements());
}

TEST_CASE("VectorTest-testProduct") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector3 = smallVector1.product(7.0);
    REQUIRE(140 == smallVector3.sumOfElements());
}

TEST_CASE("VectorTest-testL1NormalizeSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    smallVector1.l1Normalize();
    REQUIRE(1.0 == smallVector1.sumOfElements());
}

TEST_CASE("VectorTest-testL1NormalizeLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    largeVector1.l1Normalize();
    REQUIRE(1.0 == largeVector1.sumOfElements());
}

TEST_CASE("VectorTest-testL2NormSmall") {
    double data1[] = {2, 3, 4, 5, 6};
    Vector smallVector1 = Vector(data1, 5);
    double norm = smallVector1.l2Norm();
    REQUIRE(norm == sqrt(90));
}

TEST_CASE("VectorTest-testL2NormLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    double norm = largeVector1.l2Norm();
    REQUIRE(norm == sqrt(333833500));
}

TEST_CASE("VectorTest-cosineSimilaritySmall") {
    double data1[] = {2, 3, 4, 5, 6};
    double data2[] = {8, 7, 6, 5, 4};
    Vector smallVector1 = Vector(data1, 5);
    Vector smallVector2 = Vector(data2, 5);
    double similarity = smallVector1.cosineSimilarity(smallVector2);
    REQUIRE_THAT(0.8411910, Catch::Matchers::WithinAbs(similarity, 0.000001));
}

TEST_CASE("VectorTest-cosineSimilarityLarge") {
    auto *largeData1 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData1[i - 1] = i;
    }
    Vector largeVector1 = Vector(largeData1, 1000);
    auto *largeData2 = new double[1000];
    for (int i = 1; i <= 1000; i++) {
        largeData2[i - 1] = 1000 - i + 1;
    }
    Vector largeVector2 = Vector(largeData2, 1000);
    double similarity = largeVector1.cosineSimilarity(largeVector2);
    REQUIRE_THAT(0.5007497, Catch::Matchers::WithinAbs(similarity, 0.000001));
}