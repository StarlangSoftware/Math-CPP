//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#include <string>
#include <vector>
#include <cmath>
#include "DiscreteDistribution.h"

/**
 * A constructor of {@link DiscreteDistribution} class which calls its super class.
 */
DiscreteDistribution::DiscreteDistribution() : map() {
}

/**
 * The addItem method takes a String item as an input and if this map contains a mapping for the item it puts the item
 * with given value + 1, else it puts item with value of 1.
 *
 * @param item string input.
 */
void DiscreteDistribution::addItem(string item) {
    if (containsItem(item)) {
        emplace(item, find(item)->second + 1);
    } else {
        emplace(item, 1);
    }
    sum++;
}

/**
 * The removeItem method takes a String item as an input and if this map contains a mapping for the item it puts the item
 * with given value - 1, and if its value is 0, it removes the item.
 *
 * @param item String input.
 */
void DiscreteDistribution::removeItem(string item) {
    if (containsItem(item)) {
        emplace(item, find(item)->second - 1);
        if (find(item)->second == 0) {
            erase(item);
        }
    }
}

/**
 * The addDistribution method takes a {@link DiscreteDistribution} as an input and loops through the entries in this distribution
 * and if this map contains a mapping for the entry it puts the entry with its value + entry, else it puts entry with its value.
 * It also accumulates the values of entries and assigns to the sum variable.
 *
 * @param distribution {@link DiscreteDistribution} type input.
 */
void DiscreteDistribution::addDistribution(DiscreteDistribution distribution) {
    for (auto it = begin(); it != end(); it++) {
        if (containsItem(it->first)) {
            emplace(it->first, it->second + find(it->first)->second);
        } else {
            emplace(it->first, it->second);
        }
        sum += it->second;
    }
}

/**
 * The removeDistribution method takes a {@link DiscreteDistribution} as an input and loops through the entries in this distribution
 * and if this map contains a mapping for the entry it puts the entry with its key - value, else it removes the entry.
 * It also decrements the value of entry from sum and assigns to the sum variable.
 *
 * @param distribution {@link DiscreteDistribution} type input.
 */
void DiscreteDistribution::removeDistribution(DiscreteDistribution distribution) {
    for (auto it = begin(); it != end(); it++) {
        if (find(it->first)->second - it->second != 0) {
            emplace(it->first, it->second - find(it->first)->second);
        } else {
            erase(it->first);
        }
        sum -= it->second;
    }
}

/**
 * The getter for sum variable.
 *
 * @return sum.
 */
double DiscreteDistribution::getSum() {
    return sum;
}

/**
 * The containsItem method takes an item as an input and returns <tt>true</tt> if this map contains a mapping for the
 * given item.
 *
 * @param item to check.
 * @return <tt>true</tt> if this map contains a mapping for the given item.
 */
bool DiscreteDistribution::containsItem(string item) {
    return find(item) != end();
}

/**
 * The getCount method takes an item as an input returns the value to which the specified item is mapped, or {@code null}
 * if this map contains no mapping for the key.
 *
 * @param item is used to search for value.
 * @return the value to which the specified item is mapped
 */
int DiscreteDistribution::getCount(string item) {
    return find(item)->second;
}

/**
 * The getMaxItem method loops through the entries and gets the entry with maximum value.
 *
 * @return the entry with maximum value.
 */
string DiscreteDistribution::getMaxItem() {
    int max = -1;
    string maxItem = "";
    for (auto it = begin(); it != end(); it++) {
        if (it->second > max) {
            max = it->second;
            maxItem = it->first;
        }
    }
    return maxItem;
}

/**
 * Another getMaxItem method which takes an {@link vector} of Strings. It loops through the items in this {@link vector}
 * and gets the item with maximum value.
 *
 * @param includeTheseOnly {@link vector} of Strings.
 * @return the item with maximum value.
 */
string DiscreteDistribution::getMaxItem(vector<string> includeTheseOnly) {
    int max = -1;
    string maxItem = "";
    for (string item : includeTheseOnly) {
        int frequency = 0;
        if (containsItem(item)) {
            frequency = find(item)->second;
        }
        if (frequency > max) {
            max = frequency;
            maxItem = item;
        }
    }
    return maxItem;
}

/**
 * The getProbability method takes an item as an input returns the value to which the specified item is mapped over sum,
 * or 0.0 if this map contains no mapping for the key.
 *
 * @param item is used to search for probability.
 * @return the probability to which the specified item is mapped.
 */
double DiscreteDistribution::getProbability(string item) {
    if (containsItem(item)) {
        return find(item)->second / sum;
    } else {
        return 0.0;
    }
}

/**
 * The getProbabilityLaplaceSmoothing method takes an item as an input returns the smoothed value to which the specified
 * item is mapped over sum, or 1.0 over sum if this map contains no mapping for the key.
 *
 * @param item is used to search for probability.
 * @return the smoothed probability to which the specified item is mapped.
 */
double DiscreteDistribution::getProbabilityLaplaceSmoothing(string item) {
    if (containsItem(item)) {
        return (find(item)->second + 1) / (sum + size() + 1);
    } else {
        return 1.0 / (sum + size() + 1);
    }
}

/**
 * The entropy method loops through the values and calculates the entropy of these values.
 *
 * @return entropy value.
 */
double DiscreteDistribution::entropy() {
    double total = 0.0, probability;
    for (auto it = begin(); it != end(); it++) {
        probability = it->second / sum;
        total += -probability * (log(probability) / log(2));
    }
    return total;
}
