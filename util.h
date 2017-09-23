#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <attributetrajectory.h>
#include "vfkm/Grid.h"

class Util
{
public:    
  static void loadFile(std::string filename, std::vector<AttributeTrajectory>&, Grid*& g, int &noAttributes);
  static void loadFile(std::string filename, std::vector<AttributeTrajectory>&, double& xLeft, double& xRight, double& yBottom, double& yTop, double& tMin, double& tMax, int &noAttributes);
};

#endif // UTIL_H
