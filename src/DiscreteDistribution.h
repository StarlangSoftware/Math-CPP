//
// Created by Olcay Taner Yıldız on 30.11.2018.
//

#ifndef MATH_DISCRETEDISTRIBUTION_H
#define MATH_DISCRETEDISTRIBUTION_H
#include <map>
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
    double getSum();
    int getCount(const string& item);
    int getIndex(const string& item);
    vector<string> getItems();
    string getMaxItem();
    string getMaxItem(const vector<string>& includeTheseOnly);
    double getProbability(const string& item);
    map<string, double> getProbabilityDistribution();
    double getProbabilityLaplaceSmoothing(const string& item);
    double entropy();
    void serialize(ostream& outputFile);
};


#endif //MATH_DISCRETEDISTRIBUTION_H
