#ifndef ATTRIBUTETRAJECTORY_H
#define ATTRIBUTETRAJECTORY_H

#include <vector>
#include "vfkm/Vector2D.h"
#include "vfkm/Vector.h"

class AttributeTrajectory
{
private:
    std::vector<std::pair<Vector2D, VECTOR_TYPE> > points; //space x time
    std::vector<Vector>                      attributeVectors;
    int numAttributes;
public:
    AttributeTrajectory();
    AttributeTrajectory(const AttributeTrajectory &cpy);
    AttributeTrajectory(std::vector<std::pair<Vector2D,VECTOR_TYPE> > points,
                        std::vector<std::vector<VECTOR_TYPE> > attribs,
                        int numAttribs);

    size_t numberOfPoints() const { return points.size(); }

    inline std::pair<Vector2D, float> getPoint(unsigned index) const {
        return points[index];
    }
    //
    inline void addPoint(unsigned index, std::pair<Vector2D, float> newPoint, const Vector &attribs){
      points.insert((points.begin() + index), newPoint );
      attributeVectors.insert((attributeVectors.begin() + index), attribs);
    }
    //
    inline void removePoint(unsigned index) {
      points.erase(points.begin() + index);
      attributeVectors.erase(attributeVectors.begin() + index);
    }
    //
    const VECTOR_TYPE getAttribute(unsigned pointIndex, unsigned attributeIndex) const{
        Vector vec = attributeVectors.at(pointIndex);
        return vec[attributeIndex];
    }
    //
    const Vector &getAttributeVector(unsigned pointIndex) const{
        return attributeVectors.at(pointIndex);
    }
    inline int getNumAttributes() const{
        return this->numAttributes;
    }
    //
    std::string toString();
    float length();        
    //
    void getTimeRange(VECTOR_TYPE& startTime, VECTOR_TYPE& endTime);
    void   getPointAtTime(VECTOR_TYPE t, Vector2D& spatialComponent, Vector& attributeComponent);
};

void testAttributeTrajectory();

#endif // ATTRIBUTETRAJECTORY_H
