//
// Created by Olcay Taner YILDIZ on 22.12.2020.
//

#include "catch.hpp"
#include "../DiscreteDistribution.h"

TEST_CASE("DiscreteDistributionTest-testAddItem1") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE(3 == smallDistribution.getCount("item1"));
    REQUIRE(2 == smallDistribution.getCount("item2"));
    REQUIRE(1 == smallDistribution.getCount("item3"));
}

TEST_CASE("DiscreteDistributionTest-testAddItem2") {
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    for (int i = 0; i < 1000; i++){
        discreteDistribution.addItem(to_string(random() % 1000));
    }
    int count = 0;
    for (int i = 0; i < 1000; i++){
        if (discreteDistribution.containsItem(to_string(i))){
            count += discreteDistribution.getCount(to_string(i));
        }
    }
    REQUIRE(1000 == count);
}

TEST_CASE("DiscreteDistributionTest-testAddItem3") {
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    for (int i = 0; i < 1000; i++){
        discreteDistribution.addItem(to_string(random() % 1000));
    }
    for (int i = 0; i < 1000000; i++){
        discreteDistribution.addItem(to_string(random() % 1000000));
    }
    REQUIRE_THAT(discreteDistribution.size() / 1000000.0, Catch::Matchers::WithinAbs(0.632, 0.001));
}

TEST_CASE("DiscreteDistributionTest-testRemoveItem") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    smallDistribution.removeItem("item1");
    smallDistribution.removeItem("item2");
    smallDistribution.removeItem("item3");
    REQUIRE(2 == smallDistribution.getCount("item1"));
    REQUIRE(1 == smallDistribution.getCount("item2"));
}

TEST_CASE("DiscreteDistributionTest-testAddDistribution1") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    discreteDistribution.addItem("item4");
    discreteDistribution.addItem("item5");
    discreteDistribution.addItem("item5");
    discreteDistribution.addItem("item2");
    smallDistribution.addDistribution(discreteDistribution);
    REQUIRE(3 == smallDistribution.getCount("item1"));
    REQUIRE(3 == smallDistribution.getCount("item2"));
    REQUIRE(1 == smallDistribution.getCount("item3"));
    REQUIRE(1 == smallDistribution.getCount("item4"));
    REQUIRE(2 == smallDistribution.getCount("item5"));
}

TEST_CASE("DiscreteDistributionTest-testAddDistribution2") {
    DiscreteDistribution discreteDistribution1 = DiscreteDistribution();
    for (int i = 0; i < 1000; i++){
        discreteDistribution1.addItem(to_string(i));
    }
    DiscreteDistribution discreteDistribution2 = DiscreteDistribution();
    for (int i = 500; i < 1000; i++){
        discreteDistribution2.addItem(to_string(1000 + i));
    }
    discreteDistribution1.addDistribution(discreteDistribution2);
    REQUIRE(1500 == discreteDistribution1.size());
}

TEST_CASE("DiscreteDistributionTest-testRemoveDistribution") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    discreteDistribution.addItem("item1");
    discreteDistribution.addItem("item1");
    discreteDistribution.addItem("item2");
    smallDistribution.removeDistribution(discreteDistribution);
    REQUIRE(1 == smallDistribution.getCount("item1"));
    REQUIRE(1 == smallDistribution.getCount("item2"));
    REQUIRE(1 == smallDistribution.getCount("item3"));
}

TEST_CASE("DiscreteDistributionTest-testGetSum1") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE(6 == smallDistribution.getSum());
}

TEST_CASE("DiscreteDistributionTest-testGetSum2") {
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    for (int i = 0; i < 1000; i++){
        discreteDistribution.addItem(to_string(random() % 1000));
    }
    REQUIRE(1000 == discreteDistribution.getSum());
}

TEST_CASE("DiscreteDistributionTest-testGetIndex") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE(0 == smallDistribution.getIndex("item1"));
    REQUIRE(1 == smallDistribution.getIndex("item2"));
    REQUIRE(2 == smallDistribution.getIndex("item3"));
}

TEST_CASE("DiscreteDistributionTest-testContainsItem") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE_FALSE(!smallDistribution.containsItem("item1"));
    REQUIRE_FALSE(smallDistribution.containsItem("item4"));
}

TEST_CASE("DiscreteDistributionTest-testGetCount") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE(3 == smallDistribution.getCount("item1"));
    REQUIRE(2 == smallDistribution.getCount("item2"));
    REQUIRE(1 == smallDistribution.getCount("item3"));
}

TEST_CASE("DiscreteDistributionTest-testGetMaxItem1") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE("item1" == smallDistribution.getMaxItem());
}

TEST_CASE("DiscreteDistributionTest-testGetMaxItem2") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    vector<string> include;
    include.emplace_back("item2");
    include.emplace_back("item3");
    REQUIRE("item2" == smallDistribution.getMaxItem(include));
}

TEST_CASE("DiscreteDistributionTest-testGetProbability1") {
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    for (int i = 0; i < 1000; i++){
        discreteDistribution.addItem(to_string(i));
    }
    REQUIRE(0.001 == discreteDistribution.getProbability(to_string(random() % 1000)));
}

TEST_CASE("DiscreteDistributionTest-testGetProbability2") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE(0.5 == smallDistribution.getProbability("item1"));
    REQUIRE_THAT(0.333333, Catch::Matchers::WithinAbs(smallDistribution.getProbability("item2"), 0.0001));
    REQUIRE_THAT(0.166667, Catch::Matchers::WithinAbs(smallDistribution.getProbability("item3"), 0.0001));
}

TEST_CASE("DiscreteDistributionTest-getProbabilityLaplaceSmoothing1") {
    DiscreteDistribution discreteDistribution = DiscreteDistribution();
    for (int i = 0; i < 1000; i++){
        discreteDistribution.addItem(to_string(i));
    }
    REQUIRE(2.0 / 2001 == discreteDistribution.getProbabilityLaplaceSmoothing(to_string(random() % 1000)));
    REQUIRE(1.0 / 2001 == discreteDistribution.getProbabilityLaplaceSmoothing("item0"));
}

TEST_CASE("DiscreteDistributionTest-getProbabilityLaplaceSmoothing2") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE(0.4 == smallDistribution.getProbabilityLaplaceSmoothing("item1"));
    REQUIRE(0.3 == smallDistribution.getProbabilityLaplaceSmoothing("item2"));
    REQUIRE(0.2 == smallDistribution.getProbabilityLaplaceSmoothing("item3"));
    REQUIRE(0.1 == smallDistribution.getProbabilityLaplaceSmoothing("item4"));
}

TEST_CASE("DiscreteDistributionTest-testEntropy") {
    DiscreteDistribution smallDistribution = DiscreteDistribution();
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item3");
    smallDistribution.addItem("item1");
    smallDistribution.addItem("item2");
    smallDistribution.addItem("item1");
    REQUIRE_THAT(1.4591, Catch::Matchers::WithinAbs(smallDistribution.entropy(), 0.0001));
}