#include "util.h"
#include <cstdlib>
#include <cassert>
#include <fstream>
#include "vfkm/Vector2D.h"

using namespace std;

//TODO: we are assuming all trajectories are completely contained in the grid bounding box
//      still need to implement the tesselation
void Util::loadFile(std::string filename, std::vector<AttributeTrajectory>& trajectories, Grid*& grid, int &noAttributes)
{
    //clear previous data
    trajectories.clear();
    if(grid != NULL)
        delete grid;

    //
    ifstream file (filename.c_str());
    if (!file.is_open()){
      cout << "Could not open file" << endl;
      exit(0);
    }

    //parse grid info
    double xLeft, xRight, yBottom, yTop, tMin, tMax;
    file >> xLeft;
    file >> xRight;
    file >> yBottom;
    file >> yTop;
    file >> tMin;
    file >> tMax;
    
    //cout << "Grid: " << xLeft << " " << xRight << " " << yBottom << " " << yTop << " " << xResolution << " " << yResolution << endl;
    
    grid = new Grid(xLeft, yBottom, xRight - xLeft, yTop - yBottom, 10, 10);

    //line format: noTrajectories noAttributes
    int noTrajectories;
    file >> noTrajectories;
    file >> noAttributes;

    //read trajectories
    for(int trajIndex = 0 ; trajIndex < noTrajectories ; ++trajIndex){
      //parse no points in trajectory
      int noOfPoints;
      file >> noOfPoints;
      
      std::vector<std::vector<VECTOR_TYPE> > attributes;
      vector<std::pair<Vector2D,VECTOR_TYPE> > points;
      for(int i = 0 ; i < noOfPoints ; ++i){
	//
	VECTOR_TYPE x;
	VECTOR_TYPE y;
	VECTOR_TYPE t;
	file >> x >> y >> t;
	
	assert(grid->contains(x,y));
	points.push_back(make_pair(Vector2D(x,y),t));            

	//
	vector<VECTOR_TYPE> pointAttributes;
	for(int att = 0 ; att < noAttributes; ++att){
	  VECTOR_TYPE value;
	  file >> value;
	  pointAttributes.push_back(value);
	}
	attributes.push_back(pointAttributes);
      }

      //create traj
      AttributeTrajectory traj(points, attributes, noAttributes);
      trajectories.push_back(traj);
    }
}

void Util::loadFile(std::string filename, std::vector<AttributeTrajectory>& trajectories, double& xLeft, double& xRight, double& yBottom, double& yTop, double& tMin, double& tMax, int &noAttributes){
  //clear previous data
  trajectories.clear();

  //
  ifstream file (filename.c_str());
  if (!file.is_open()){
    cout << "Could not open file" << endl;
    exit(0);
  }

  //parse grid info
  file >> xLeft;
  file >> xRight;
  file >> yBottom;
  file >> yTop;
  file >> tMin;
  file >> tMax;
    
  //line format: noTrajectories noAttributes
  int noTrajectories;
  file >> noTrajectories;
  file >> noAttributes;

  //read trajectories
  for(int trajIndex = 0 ; trajIndex < noTrajectories ; ++trajIndex){
    //parse no points in trajectory
    int noOfPoints;
    file >> noOfPoints;
      
    std::vector<std::vector<VECTOR_TYPE> > attributes;
    vector<std::pair<Vector2D,VECTOR_TYPE> > points;
    for(int i = 0 ; i < noOfPoints ; ++i){
      //
      VECTOR_TYPE x;
      VECTOR_TYPE y;
      VECTOR_TYPE t;
      file >> x >> y >> t;

      //assuming here that the point is inside the bounding box
      points.push_back(make_pair(Vector2D(x,y),t));            

      //
      vector<VECTOR_TYPE> pointAttributes;
      for(int att = 0 ; att < noAttributes; ++att){
	VECTOR_TYPE value;
	file >> value;
	pointAttributes.push_back(value);
      }
      attributes.push_back(pointAttributes);
    }

    //create traj
    AttributeTrajectory traj(points, attributes, noAttributes);
    trajectories.push_back(traj);
  }
}
