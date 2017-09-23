#include "attributetrajectory.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "util.h"
#include "vfkm/Optimizer.h"
//#include "randomsampler.h"

using namespace std;

//string initializationToString(AFKMInitializationMethod x){
//    if(x == ALL_RANDOM)
//        return "ALL_RANDOM";
//    else if(x == KMEANS_PLUS_PLUS)
//        return "KMEANS_PLUS_PLUS";
//    else if(x == FIRST_RANDOM_WORST_MATCH)
//        return "FIRST_RANDOM_WORST_MATCH";
//    else if(x == FROM_FILE)
//        return "FROM_FILE";
//    else
//        assert(false);
//}

//void sampleExperiment(string directory, string name, time_t randomSeed,AFKMInitializationMethod initializationMethod){
//    string filename = directory + "/" + name;
//
//    //
//    srand(randomSeed);
//    Grid* grid = NULL;
//    std::vector<AttributeTrajectory> trajectories;
//    int numberOfAttributes;
//    Util::loadFile(filename,trajectories, grid, numberOfAttributes);
//    int numberOfCurves = trajectories.size();
//    cout << "Loaded " << trajectories.size() << " trajectories.";
//    cout << "    Number of Attributes " << numberOfAttributes;
//
//    //
//    vector<vector<Vector> > attributeFields;
//    vector<double> errorPerAttribute;
//    int numberOfFields = 4;
//    unsigned short mapCurveToField[numberOfCurves];
//    float mapCurveToError[numberOfCurves];
//    float smoothnessWeight = 0.05;
//
//    //
//    double totalError = 0.0;
//    int numberOfIterations = 0;
//    Optimizer::optimizeAFKM(*grid, numberOfFields,numberOfAttributes,trajectories, attributeFields, mapCurveToField, mapCurveToError, totalError, numberOfIterations,errorPerAttribute, smoothnessWeight,100,initializationMethod);
//
//    //
//    cout << "Result: Total Error " << totalError << endl;
//    vector<int> countPerField(numberOfFields,0);
//    QFile resultFile(directory + "/cluster_result_synthetic.txt");
//    QTextStream stream(&resultFile);
//    if (!resultFile.open(QIODevice::WriteOnly)){
//        cout << "Could not write file" << endl;
//        assert(false);
//    }
//
//    for(int i = 0 ; i < numberOfCurves ; ++i){
//        stream << mapCurveToField[i] << "\n";
//        countPerField[mapCurveToField[i]] += 1;
//    }
//
//    for(int i = 0 ; i < numberOfFields ; ++i){
//        cout << "    Field " << i << " has " << countPerField[i] << " trajectories" << endl;
//    }
//
//    for(int d = 0 ; d < numberOfAttributes ; ++d){
//        cout << " Attribute " << d << " error = " << errorPerAttribute.at(d) << endl;
//    }
//
//    //cout << endl << "Fields:" << endl;
//    for(int i = 0 ; i < numberOfFields ; ++i){
//        //cout << "   Field " << i << ":" << endl;
//        vector<Vector>& currentField = attributeFields.at(i);
//        for(int d = 0 ; d < numberOfAttributes ; ++d){
//            //cout << "      Dimension " << d << ":" << endl;
//            QString filename;
//            QTextStream filenameStream(&filename);
//            filenameStream << directory << "/field_" << i << "_" << d << ".txt";
//            QFile fieldFile(filename);
//            if (!fieldFile.open(QIODevice::WriteOnly)) {
//                cout << "Could not open file" << endl;
//                return;
//            }
//            QTextStream stream(&fieldFile);
//
//            //
//            Vector& dimensionField = currentField.at(d);
//#if 0
//            int numVertices = dimensionField.getDimension();
//            for(int v = 0 ; v < numVertices ; ++v){
//                cout << "         Value " << v << " = " << dimensionField[v] << endl;
//            }
//#elif 0
//            int resX = grid->getResolutionX();
//            int resY = grid->getResolutionY();
//            QString str;
//            QTextStream stream(&str);
//            stream << "[";
//            for(int i = 0 ; i < resX ; ++i){
//                stream << "[";
//                for(int j = 0 ; j < resY ; ++j){
//                    stream << dimensionField[i * resY + j];
//                    if(j < (resY - 1))
//                        stream << ",";
//                }
//                stream << "]";
//            }
//            stream << "]";
//            cout << str.toStdString() << endl;
//#else
//            int resX = grid->getResolutionX();
//            int resY = grid->getResolutionY();
//            for(int i = 0 ; i < resX ; ++i){
//                for(int j = 0 ; j < resY ; ++j){
//                    stream << dimensionField[i * resY + j];
//                    if(j < (resY - 1))
//                        stream << ",";
//                }
//                stream << "\n";
//            }
//#endif
//        }
//    }
//
//}

void testLoadFile(){
  std::vector<AttributeTrajectory> trajectories;
  Grid* g = NULL;
  string filename("../AFKM/data/atlantic.txt");
  int noAttributes;
  Util::loadFile(filename, trajectories, g, noAttributes);

  cout << "Loaded " << trajectories.size() << " trajectories" << endl;
  int numTrajectories = trajectories.size();
  for(int i = 0 ; i < numTrajectories ; ++i){
    cout << trajectories.at(i).toString() << endl;
  }
}


void initExperiment(string filename, Grid*& g, std::vector<AttributeTrajectory>& curves,
		    int gridXResolution, int gridYResolution, int& noAttributes){
  double xMin;
  double xMax;
  double yMin;
  double yMax;
  double tMin;
  double tMax;
  Util::loadFile(filename, curves, xMin, xMax, yMin, yMax, tMin, tMax, noAttributes);

  //
  g = new Grid(xMin,yMin,xMax-xMin,yMax-yMin,gridXResolution, gridYResolution);
}

int main(int argc, char *argv[])
{
  int rightNumberOfParameters = 7;
  
  if(argc != rightNumberOfParameters){
    //print usage
    cout << "./afkm trajectoryFile gridXResolution gridYResolution numberOfClusterFields smoothnessWeight outputDirectory" << endl;
    return 0;
  }

  //set parameters
  string filename(argv[1]);
  int    gridXResolution = atoi(argv[2]);
  int    gridYResolution = atoi(argv[3]);
  int    numberOfFields = atoi(argv[4]);
  float  smoothnessWeight = atof(argv[5]);
  string outputDirectory(argv[6]);

  //load files
  cout << "Loading Files..." << endl;
  std::vector<AttributeTrajectory> curves;
  Grid* grid = NULL;
  int numberOfAttributes;
  initExperiment(filename, grid, curves, gridXResolution, gridYResolution, numberOfAttributes);

  int numberOfCurves = curves.size();
  cout << "Loaded " << numberOfCurves << " trajectories" << endl;

  //optimize
  cout << "Optimizing..." << endl;
  vector<vector<Vector> > attributeFields;
  vector<double> errorPerAttribute;
  unsigned short mapCurveToField[numberOfCurves];
  float mapCurveToError[numberOfCurves];
  
  //
  double totalError = 0.0;
  int numberOfIterations = 0;
  AFKMInitializationMethod initializationMethod = KMEANS_PLUS_PLUS;
  Optimizer::optimizeAFKM(*grid, numberOfFields,numberOfAttributes, curves, attributeFields, mapCurveToField, mapCurveToError, totalError, numberOfIterations,errorPerAttribute, smoothnessWeight,100,initializationMethod);

  //
  cout << "Result: Total Error " << totalError << endl;
  vector<int> countPerField(numberOfFields,0);
  stringstream ss;
  ss << outputDirectory + "/cluster_result_synthetic.txt";
  ofstream resultFile(ss.str());
  if (!resultFile.is_open()){
    cout << "Could not write file" << endl;
    assert(false);
  }
  
  for(int i = 0 ; i < numberOfCurves ; ++i){
    resultFile << mapCurveToField[i] << "\n";
    countPerField[mapCurveToField[i]] += 1;
  }

  for(int i = 0 ; i < numberOfFields ; ++i){
    cout << "    Field " << i << " has " << countPerField[i] << " trajectories" << endl;
  }

  for(int d = 0 ; d < numberOfAttributes ; ++d){
    cout << " Attribute " << d << " error = " << errorPerAttribute.at(d) << endl;
  }
  //cout << endl << "Fields:" << endl;
  for(int i = 0 ; i < numberOfFields ; ++i){
    //cout << "   Field " << i << ":" << endl;
    vector<Vector>& currentField = attributeFields.at(i);
    for(int d = 0 ; d < numberOfAttributes ; ++d){
      //cout << "      Dimension " << d << ":" << endl;
      stringstream filename;
      filename << outputDirectory << "/field_" << i << "_" << d << ".txt";
      ofstream fieldFile(filename.str().c_str());
      if (!fieldFile.is_open()) {
	cout << "Could not open file result file: " << filename.str() << endl;
	assert(false);
      }

      //
      Vector& dimensionField = currentField.at(d);

      int resX = grid->getResolutionX();
      int resY = grid->getResolutionY();
      for(int i = 0 ; i < resX ; ++i){
	for(int j = 0 ; j < resY ; ++j){
	  fieldFile << dimensionField[i * resY + j];
	  if(j < (resY - 1))
	    fieldFile << ",";
	}
	fieldFile << "\n";
      }
    }
  }

  //
  return 0;
}
