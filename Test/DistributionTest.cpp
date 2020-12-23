//
// Created by Olcay Taner YILDIZ on 22.12.2020.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../Distribution.h"

TEST_CASE("DistributionTest-testZNormal") {
    REQUIRE_THAT(0.5, Catch::Matchers::WithinAbs(Distribution::zNormal(0.0), 0.0));
    REQUIRE_THAT(0.69146, Catch::Matchers::WithinAbs(Distribution::zNormal(0.5), 0.00001));
    REQUIRE_THAT(0.84134, Catch::Matchers::WithinAbs(Distribution::zNormal(1.0), 0.00001));
    REQUIRE_THAT(0.93319, Catch::Matchers::WithinAbs(Distribution::zNormal(1.5), 0.00001));
    REQUIRE_THAT(0.97725, Catch::Matchers::WithinAbs(Distribution::zNormal(2.0), 0.00001));
    REQUIRE_THAT(0.99379, Catch::Matchers::WithinAbs(Distribution::zNormal(2.5), 0.00001));
    REQUIRE_THAT(0.99865, Catch::Matchers::WithinAbs(Distribution::zNormal(3.0), 0.00001));
    REQUIRE_THAT(0.99977, Catch::Matchers::WithinAbs(Distribution::zNormal(3.5), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(0.5), Catch::Matchers::WithinAbs(Distribution::zNormal(-0.5), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(1.0), Catch::Matchers::WithinAbs(Distribution::zNormal(-1.0), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(1.5), Catch::Matchers::WithinAbs(Distribution::zNormal(-1.5), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(2.0), Catch::Matchers::WithinAbs(Distribution::zNormal(-2.0), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(2.5), Catch::Matchers::WithinAbs(Distribution::zNormal(-2.5), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(3.0), Catch::Matchers::WithinAbs(Distribution::zNormal(-3.0), 0.00001));
    REQUIRE_THAT(1 - Distribution::zNormal(3.5), Catch::Matchers::WithinAbs(Distribution::zNormal(-3.5), 0.00001));
}

TEST_CASE("DistributionTest-testZInverse") {
    REQUIRE_THAT(0.0, Catch::Matchers::WithinAbs(Distribution::zInverse(0.5), 0.00001));
    REQUIRE_THAT(0.841621, Catch::Matchers::WithinAbs(Distribution::zInverse(0.8), 0.00001));
    REQUIRE_THAT(1.281552, Catch::Matchers::WithinAbs(Distribution::zInverse(0.9), 0.00001));
    REQUIRE_THAT(1.644854, Catch::Matchers::WithinAbs(Distribution::zInverse(0.95), 0.00001));
    REQUIRE_THAT(2.053749, Catch::Matchers::WithinAbs(Distribution::zInverse(0.98), 0.00001));
    REQUIRE_THAT(2.326348, Catch::Matchers::WithinAbs(Distribution::zInverse(0.99), 0.00001));
    REQUIRE_THAT(2.575829, Catch::Matchers::WithinAbs(Distribution::zInverse(0.995), 0.00001));
    REQUIRE_THAT(2.878162, Catch::Matchers::WithinAbs(Distribution::zInverse(0.998), 0.00001));
    REQUIRE_THAT(3.090232, Catch::Matchers::WithinAbs(Distribution::zInverse(0.999), 0.00001));
}

TEST_CASE("DistributionTest-testChiSquare") {
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::chiSquare(3.841, 1), 0.0001));
    REQUIRE_THAT(0.005, Catch::Matchers::WithinAbs(Distribution::chiSquare(7.879, 1), 0.0001));
    REQUIRE_THAT(0.95, Catch::Matchers::WithinAbs(Distribution::chiSquare(3.940, 10), 0.0001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::chiSquare(18.307, 10), 0.0001));
    REQUIRE_THAT(0.995, Catch::Matchers::WithinAbs(Distribution::chiSquare(2.156, 10), 0.0001));
    REQUIRE_THAT(0.005, Catch::Matchers::WithinAbs(Distribution::chiSquare(25.188, 10), 0.0001));
    REQUIRE_THAT(0.95, Catch::Matchers::WithinAbs(Distribution::chiSquare(77.929, 100), 0.0001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::chiSquare(124.342, 100), 0.0001));
    REQUIRE_THAT(0.995, Catch::Matchers::WithinAbs(Distribution::chiSquare(67.328, 100), 0.0001));
    REQUIRE_THAT(0.005, Catch::Matchers::WithinAbs(Distribution::chiSquare(140.169, 100), 0.0001));
}

TEST_CASE("DistributionTest-testChiSquareInverse") {
    REQUIRE_THAT(2.706, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.1, 1), 0.001));
    REQUIRE_THAT(6.635, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.01, 1), 0.001));
    REQUIRE_THAT(4.865, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.9, 10), 0.001));
    REQUIRE_THAT(15.987, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.1, 10), 0.001));
    REQUIRE_THAT(2.558, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.99, 10), 0.001));
    REQUIRE_THAT(23.209, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.01, 10), 0.001));
    REQUIRE_THAT(82.358, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.9, 100), 0.001));
    REQUIRE_THAT(118.498, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.1, 100), 0.001));
    REQUIRE_THAT(70.065, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.99, 100), 0.001));
    REQUIRE_THAT(135.807, Catch::Matchers::WithinAbs(Distribution::chiSquareInverse(0.01, 100), 0.001));
}

TEST_CASE("DistributionTest-testFDistribution") {
    REQUIRE_THAT(0.1, Catch::Matchers::WithinAbs(Distribution::fDistribution(39.86346, 1, 1), 0.00001));
    REQUIRE_THAT(0.1, Catch::Matchers::WithinAbs(Distribution::fDistribution(2.32260, 10, 10), 0.00001));
    REQUIRE_THAT(0.1, Catch::Matchers::WithinAbs(Distribution::fDistribution(1.79384, 20, 20), 0.00001));
    REQUIRE_THAT(0.1, Catch::Matchers::WithinAbs(Distribution::fDistribution(1.60648, 30, 30), 0.00001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::fDistribution(161.4476, 1, 1), 0.00001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::fDistribution(2.9782, 10, 10), 0.00001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::fDistribution(2.1242, 20, 20), 0.00001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::fDistribution(1.8409, 30, 30), 0.00001));
    REQUIRE_THAT(0.01, Catch::Matchers::WithinAbs(Distribution::fDistribution(4052.181, 1, 1), 0.00001));
    REQUIRE_THAT(0.01, Catch::Matchers::WithinAbs(Distribution::fDistribution(4.849, 10, 10), 0.00001));
    REQUIRE_THAT(0.01, Catch::Matchers::WithinAbs(Distribution::fDistribution(2.938, 20, 20), 0.00001));
    REQUIRE_THAT(0.01, Catch::Matchers::WithinAbs(Distribution::fDistribution(2.386, 30, 30), 0.00001));
}

TEST_CASE("DistributionTest-testFDistributionInverse") {
    REQUIRE_THAT(3.818, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.01, 5, 26), 0.001));
    REQUIRE_THAT(15.1010, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.025, 4, 3), 0.001));
    REQUIRE_THAT(2.19535, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.1, 8, 13), 0.001));
    REQUIRE_THAT(2.29871, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.1, 3, 27), 0.001));
    REQUIRE_THAT(3.4381, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.05, 8, 8), 0.001));
    REQUIRE_THAT(2.6283, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.05, 6, 19), 0.001));
    REQUIRE_THAT(3.3120, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.025, 9, 13), 0.001));
    REQUIRE_THAT(3.7505, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.025, 3, 23), 0.001));
    REQUIRE_THAT(4.155, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.01, 12, 12), 0.001));
    REQUIRE_THAT(6.851, Catch::Matchers::WithinAbs(Distribution::fDistributionInverse(0.01, 1, 120), 0.001));
}

TEST_CASE("DistributionTest-testTDistribution") {
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::tDistribution(6.314, 1), 0.0001));
    REQUIRE_THAT(0.005, Catch::Matchers::WithinAbs(Distribution::tDistribution(63.656, 1), 0.0001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::tDistribution(1.812, 10), 0.0001));
    REQUIRE_THAT(0.01, Catch::Matchers::WithinAbs(Distribution::tDistribution(2.764, 10), 0.0001));
    REQUIRE_THAT(0.005, Catch::Matchers::WithinAbs(Distribution::tDistribution(3.169, 10), 0.0001));
    REQUIRE_THAT(0.001, Catch::Matchers::WithinAbs(Distribution::tDistribution(4.144, 10), 0.0001));
    REQUIRE_THAT(0.05, Catch::Matchers::WithinAbs(Distribution::tDistribution(1.725, 20), 0.0001));
    REQUIRE_THAT(0.01, Catch::Matchers::WithinAbs(Distribution::tDistribution(2.528, 20), 0.0001));
    REQUIRE_THAT(0.005, Catch::Matchers::WithinAbs(Distribution::tDistribution(2.845, 20), 0.0001));
    REQUIRE_THAT(0.001, Catch::Matchers::WithinAbs(Distribution::tDistribution(3.552, 20), 0.0001));
}

TEST_CASE("DistributionTest-testTDistributionInverse") {
    REQUIRE_THAT(2.947, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.005, 15), 0.001));
    REQUIRE_THAT(1.717, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.05, 22), 0.001));
    REQUIRE_THAT(3.365, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.01, 5), 0.001));
    REQUIRE_THAT(3.922, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.0005, 18), 0.001));
    REQUIRE_THAT(3.467, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.001, 24), 0.001));
    REQUIRE_THAT(6.314, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.05, 1), 0.001));
    REQUIRE_THAT(2.306, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.025, 8), 0.001));
    REQUIRE_THAT(3.646, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.001, 17), 0.001));
    REQUIRE_THAT(3.373, Catch::Matchers::WithinAbs(Distribution::tDistributionInverse(0.0005, 120), 0.001));
}