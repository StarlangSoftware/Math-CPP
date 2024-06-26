//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <ostream>
#include "DiscreteDistribution.h"

/**
 * A constructor of DiscreteDistribution class which calls its super class.
 */
DiscreteDistribution::DiscreteDistribution() : map() {
}

/**
 * The addItem method takes a String item as an input and if this map contains a mapping for the item it puts the item
 * with given value + 1, else it puts item with value of 1.
 *
 * @param item string input.
 */
void DiscreteDistribution::addItem(const string& item) {
    if (contains(item)) {
        insert_or_assign(item, find(item)->second + 1);
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
void DiscreteDistribution::removeItem(const string& item) {
    if (contains(item)) {
        insert_or_assign(item, find(item)->second - 1);
        if (find(item)->second == 0) {
            erase(item);
        }
        sum--;
    }
}

/**
 * The addDistribution method takes a DiscreteDistribution as an input and loops through the entries in this distribution
 * and if this map contains a mapping for the entry it puts the entry with its value + entry, else it puts entry with its value.
 * It also accumulates the values of entries and assigns to the sum variable.
 *
 * @param distribution DiscreteDistribution type input.
 */
void DiscreteDistribution::addDistribution(const DiscreteDistribution& distribution) {
    for (auto &it : distribution) {
        if (contains(it.first)) {
            insert_or_assign(it.first, it.second + find(it.first)->second);
        } else {
            emplace(it.first, it.second);
        }
        sum += it.second;
    }
}

/**
 * The removeDistribution method takes a DiscreteDistribution as an input and loops through the entries in this distribution
 * and if this map contains a mapping for the entry it puts the entry with its key - value, else it removes the entry.
 * It also decrements the value of entry from sum and assigns to the sum variable.
 *
 * @param distribution DiscreteDistribution type input.
 */
void DiscreteDistribution::removeDistribution(const DiscreteDistribution& distribution) {
    for (auto &it : distribution) {
        if (find(it.first)->second - it.second != 0) {
            insert_or_assign(it.first, find(it.first)->second - it.second);
        } else {
            erase(it.first);
        }
        sum -= it.second;
    }
}

/**
 * The getter for sum variable.
 *
 * @return sum.
 */
double DiscreteDistribution::getSum() const{
    return sum;
}

/**
 * The getCount method takes an item as an input returns the value to which the specified item is mapped, or {@code null}
 * if this map contains no mapping for the key.
 *
 * @param item is used to search for value.
 * @return the value to which the specified item is mapped
 */
int DiscreteDistribution::getCount(const string& item) const{
    return find(item)->second;
}

/**
 * The getMaxItem method loops through the entries and gets the entry with maximum value.
 *
 * @return the entry with maximum value.
 */
string DiscreteDistribution::getMaxItem() const{
    int max = -1;
    string maxItem;
    for (auto &it : *this) {
        if (it.second > max) {
            max = it.second;
            maxItem = it.first;
        }
    }
    return maxItem;
}

/**
 * Another getMaxItem method which takes an vector of Strings. It loops through the items in this vector
 * and gets the item with maximum value.
 *
 * @param includeTheseOnly vector of Strings.
 * @return the item with maximum value.
 */
string DiscreteDistribution::getMaxItem(const vector<string>& includeTheseOnly) const{
    int max = -1;
    string maxItem;
    for (const string &item : includeTheseOnly) {
        int frequency = 0;
        if (contains(item)) {
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
double DiscreteDistribution::getProbability(const string& item) const{
    if (contains(item)) {
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
double DiscreteDistribution::getProbabilityLaplaceSmoothing(const string& item) const{
    if (contains(item)) {
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
double DiscreteDistribution::entropy() const{
    double total = 0.0, probability;
    for (auto &it : *this) {
        probability = it.second / sum;
        total += -probability * (log(probability) / log(2));
    }
    return total;
}

/**
 * The getIndex method takes an item as an input and returns the index of given item.
 *
 * @param item to search for index.
 * @return index of given item.
 */
int DiscreteDistribution::getIndex(const string& item) const{
    int i = 0;
    for (auto &it : *this) {
        if (it.first == item){
            return i;
        }
        i++;
    }
    return -1;
}

/**
 * Returns the items in the distribution as a vector
 * @return Items in the distribution as a vector
 */
vector<string> DiscreteDistribution::getItems() const{
    vector<string> result;
    for (auto &it : *this){
        result.push_back(it.first);
    }
    return result;
}

/**
 * Write the distribution to an output file.
 * @param outputFile OUtput file
 */
void DiscreteDistribution::serialize(ostream &outputFile) {
    outputFile << sum << "\n";
    outputFile << this->size() << "\n";
    for (auto &it : *this) {
        outputFile << it.first << "\n";
        outputFile << it.second << "\n";
    }
}

/**
 * Reads the distribution from an input file
 * @param inputFile Input file
 */
DiscreteDistribution::DiscreteDistribution(ifstream &inputFile) {
    int size;
    string item;
    int count;
    sum = 0;
    inputFile >> size;
    for (int i = 0; i < size; i++){
        inputFile >> item;
        inputFile >> count;
        sum += count;
        emplace(item, count);
    }
}

/**
 * Returns the distribution as a probability distribution
 * @return Probability distribution
 */
map<string, double> DiscreteDistribution::getProbabilityDistribution() const{
    map<string, double> result;
    for (auto &it : *this) {
        result.insert_or_assign(it.first, getProbability(it.first));
    }
    return result;
}
