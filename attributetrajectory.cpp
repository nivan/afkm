#include "attributetrajectory.h"
#include <cassert>
#include <sstream>
#include <iostream>
#define DEBUG

using namespace std;

AttributeTrajectory::AttributeTrajectory():
    numAttributes(0)
{}

AttributeTrajectory::AttributeTrajectory(const AttributeTrajectory &cpy):
    numAttributes(cpy.getNumAttributes()){
    int numPoints = cpy.numberOfPoints();
    for(int i = 0 ; i < numPoints ; ++i){
        this->points.push_back(cpy.getPoint(i));
    }
    //
    for(int i = 0 ; i < numPoints ; ++i){
        Vector attribVector(numAttributes);
        const Vector& currentAttribs = cpy.getAttributeVector(i);
        attribVector.setValues(currentAttribs);
        this->attributeVectors.push_back(attribVector);
    }
}

AttributeTrajectory::AttributeTrajectory(std::vector<std::pair<Vector2D,VECTOR_TYPE> > points,
                                         std::vector<std::vector<VECTOR_TYPE> > attribs,
                                         int numAttribs):
    numAttributes(numAttribs)
{
#ifdef DEBUG
    assert(points.size() == attribs.size());
#endif
    //
    this->points.assign(points.begin(),points.end());
    int numPoints = this->points.size();
    //
    for(int i = 0 ; i < numPoints ; ++i){
        Vector attribVector(numAttributes);
        vector<VECTOR_TYPE>& currentAttribs = attribs.at(i);
#ifdef DEBUG
	int currentAttribsSize = currentAttribs.size();
        assert(currentAttribsSize == numAttribs);
#endif
        attribVector.setValues((VECTOR_TYPE*)(&currentAttribs[0]));
	attributeVectors.push_back(attribVector);
    }
}

std::string AttributeTrajectory::toString(){
    stringstream ss;

    ss << "Curve:\n";

    int numberOfPoints = points.size();

    for(int i = 0 ; i < numberOfPoints ; ++i){
        pair<Vector2D,float> currentPoint = points.at(i);
        Vector& attribVec = attributeVectors.at(i);

        ss << "   space = " << currentPoint.first.toString() << ", time = " << currentPoint.second << ", attribute = " << attribVec.toString() << ";\n";
    }

    return ss.str();
}

float AttributeTrajectory::length(){
    float total = 0;
    int numberPoints = numberOfPoints();

    for(int i = 0 ; i < numberPoints - 1 ; ++i){
        Vector2D p0 = points.at(i).first;
        Vector2D p1 = points.at(i+1).first;

        p1.subtract(p0);

        total += p1.length();
    }

    return total;
}

void AttributeTrajectory::getTimeRange(VECTOR_TYPE& startTime, VECTOR_TYPE& endTime){
    if(points.size() == 0){
        startTime = 1;
        endTime   = -1;
    }


    //
    startTime = points.front().second;
    endTime   = points.back().second;
}

void AttributeTrajectory::getPointAtTime(double t, Vector2D &spatialComponent, Vector &attributeComponent){
    double startTime = points.front().second;
    double endTime   = points.back().second;

    //
    if(t < startTime || t > endTime)
        return;

    //
    int numPoints = points.size();
    for(int i = 0 ; i < (numPoints - 1) ; ++i){
        std::pair<Vector2D, VECTOR_TYPE>& currentSample = points.at(i);
        VECTOR_TYPE currentTime = currentSample.second;

        std::pair<Vector2D, VECTOR_TYPE>&    nextSample = points.at(i+1);
        VECTOR_TYPE nextTime = nextSample.second;

        //
        if(currentTime > t || t > nextTime)
            continue;

        //found interval
        Vector2D currentPoint = currentSample.first;
        Vector2D nextPoint    = nextSample.first;

        //
        float lambda = (t - currentTime) / (nextTime - currentTime);
        currentPoint.scale(1.0 - lambda);
        nextPoint.scale(lambda);
        spatialComponent = currentPoint + nextPoint;

        //
        Vector currentAttributeVector = attributeVectors.at(i);
        Vector    nextAttributeVector = attributeVectors.at(i+1);
        attributeComponent = ((1.0 - lambda) * currentAttributeVector) + (lambda * nextAttributeVector);
    }
}

void testAttributeTrajectory(){
    vector<pair<Vector2D,VECTOR_TYPE> > points;
    //create points
    points.push_back(make_pair(Vector2D(0,0),0));
    points.push_back(make_pair(Vector2D(1,0),0.25));
    points.push_back(make_pair(Vector2D(1,1),0.5));
    points.push_back(make_pair(Vector2D(0,1),0.75));
    points.push_back(make_pair(Vector2D(0,0),1));
    //create attributes
    vector<vector<VECTOR_TYPE> >attributes;
    //
    vector<VECTOR_TYPE> attPoint1;
    attPoint1.push_back(-1);
    attributes.push_back(attPoint1);
    //
    vector<VECTOR_TYPE> attPoint2;
    attPoint2.push_back(10);
    attributes.push_back(attPoint2);
    //
    vector<VECTOR_TYPE> attPoint3;
    attPoint3.push_back(21);
    attributes.push_back(attPoint3);
    //
    vector<VECTOR_TYPE> attPoint4;
    attPoint4.push_back(81);
    attributes.push_back(attPoint4);
    //
    vector<VECTOR_TYPE> attPoint5;
    attPoint5.push_back(-10);
    attributes.push_back(attPoint5);

    //
    AttributeTrajectory traj(points,attributes,1);
    cout << traj.toString() << endl;

    //
    VECTOR_TYPE t = 0.25;
    Vector2D interpVec;
    Vector   interpAttribute;
    traj.getPointAtTime(t,interpVec,interpAttribute);
    cout << "Spatial Component = " << interpVec.toString() << endl;
    cout << "Time = " << t << endl;
    cout << "Attribute Component = " << interpAttribute.toString() << endl;
}
