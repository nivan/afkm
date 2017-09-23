#include "randomsampler.h"
#include <ctime>
#include <iostream>
#include <cstdlib>

using namespace std;

int RandomSampler::seed = time(NULL);

void RandomSampler::randomSample(int numSamples,int lowerBound,int upperBound,vector<int>& result){

}

void RandomSampler::randomSampleProbabilityDistribution(const std::vector<std::pair<int,float> >&distribution,int numSamples,std::vector<int>& result){
    result.clear();

    vector<float> cummulativeProbabilities;
    double currentProb = 0.0;
    int numElements = distribution.size();
    for(int i = 0 ; i < numElements ; ++i){
        cummulativeProbabilities.push_back(currentProb);
        currentProb += distribution.at(i).second;
        if(i == (numElements -1))
            cummulativeProbabilities.push_back(currentProb);
    }

//    for(int i = 0 ; i < cummulativeProbabilities.size() ; ++i){
//        cout << "Element " << i << " = " << cummulativeProbabilities.at(i) << endl;
//    }

    for(int i = 0 ; i < numSamples ; ++i){
        double t = (1.0*rand()) / RAND_MAX;
        auto it = lower_bound(cummulativeProbabilities.begin(), cummulativeProbabilities.end(),t);
        int index = (it - cummulativeProbabilities.begin()) - 1;
        index = max(index,0);
        //cout << "Lower Bound = " << *it << " Index " << index << endl;
        result.push_back(distribution.at(index).first);
    }
}

int RandomSampler::randomSampleProbabilityDistribution(const std::vector<std::pair<int,float> >&distribution){
    vector<float> cummulativeProbabilities;
    double currentProb = 0.0;
    int numElements = distribution.size();
    for(int i = 0 ; i < numElements ; ++i){
        cummulativeProbabilities.push_back(currentProb);
        currentProb += distribution.at(i).second;
        if(i == (numElements -1))
            cummulativeProbabilities.push_back(currentProb);
    }

    double t = (1.0*rand()) / RAND_MAX;
    auto it = lower_bound(cummulativeProbabilities.begin(), cummulativeProbabilities.end(),t);
    int index = (it - cummulativeProbabilities.begin()) - 1;
    index = max(index,0);
    return distribution.at(index).first;
}
