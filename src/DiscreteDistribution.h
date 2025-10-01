//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#ifndef MATH_DISCRETEDISTRIBUTION_H
#define MATH_DISCRETEDISTRIBUTION_H
#include <map>
#include <string>
#include <vector>
using namespace std;

class DiscreteDistribution : public map<string, int>{
private:
    double sum = 0;
public:
    DiscreteDistribution();
    explicit DiscreteDistribution(ifstream& inputFile);
    void addItem(const string& item);
    void removeItem(const string& item);
    void addDistribution(const DiscreteDistribution& distribution);
    void removeDistribution(const DiscreteDistribution& distribution);
    double getSum() const;
    int getCount(const string& item) const;
    int getIndex(const string& item) const;
    vector<string> getItems() const;
    string getMaxItem() const;
    string getMaxItem(const vector<string>& includeTheseOnly) const;
    double getProbability(const string& item) const;
    map<string, double> getProbabilityDistribution() const;
    double getProbabilityLaplaceSmoothing(const string& item) const;
    double entropy() const;
    void serialize(ostream& outputFile);
};


#endif //MATH_DISCRETEDISTRIBUTION_H
