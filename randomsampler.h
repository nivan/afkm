#ifndef RANDOMSAMPLER_H
#define RANDOMSAMPLER_H

#include <vector>

class RandomSampler
{
private:
    static int seed;

public:
    static void randomSample(int numSamples,int lowerBound,int upperBound,std::vector<int>&);

    static int  randomSampleProbabilityDistribution(const std::vector<std::pair<int,float> >&);
    static void randomSampleProbabilityDistribution(const std::vector<std::pair<int,float> >&,int numSamles,std::vector<int>&);

};

#endif // RANDOMSAMPLER_H
