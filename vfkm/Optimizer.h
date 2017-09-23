#ifndef OPTIMIZER_H
#define OPTIMIZER_H

//#include <vector>
#include "Grid.h"
#include "attributetrajectory.h"

struct ProblemSettings
{
    Grid &grid;
    const std::vector<int> &curveIndices;
    const std::vector<CurveDescription> &curve_descriptions;
    float totalCurveLength;
    float smoothnessWeight;
    ProblemSettings(Grid &g,
                    const std::vector<int> &i,
                    const std::vector<CurveDescription> &cd,
                    float tcl, float sw):
        grid(g),
        curveIndices(i),
        curve_descriptions(cd),
        totalCurveLength(tcl),
        smoothnessWeight(sw) {}
};

struct AFKMProblemSettings
{
    Grid &grid;
    const std::vector<int> &curveIndices;
    const std::vector<AFKMCurveDescription> &curve_descriptions;
    VECTOR_TYPE totalTimeSpan;
    VECTOR_TYPE smoothnessWeight;
    AFKMProblemSettings(Grid &g,
                        const std::vector<int> &i,
                        const std::vector<AFKMCurveDescription> &cd,
                        VECTOR_TYPE tcl, VECTOR_TYPE sw):
        grid(g),
        curveIndices(i),
        curve_descriptions(cd),
        totalTimeSpan(tcl),
        smoothnessWeight(sw) {}
};

enum AFKMInitializationMethod{
    ALL_RANDOM,
    FIRST_RANDOM_WORST_MATCH,
    KMEANS_PLUS_PLUS,
    FROM_FILE
};

class Optimizer
{
public:
    Optimizer(int size);
    ~Optimizer();

    static void multiplyByA(const Vector& x, Vector &resultX, Vector &diagM, ProblemSettings &prob);

    static void multiplyByA(const Vector& x, const Vector& y, Vector &resultX, Vector &resultY, Grid &grid,
                            const std::vector<int> &curveIndices,
                            std::vector< std::vector<Intersection> > &mapCurveToConstraints,
                            float totalCurveLength, float smoothnessWeight);

    static void multiplyByAWithoutWeights(Vector& x, Vector& y, Vector &resultX,
                                          Vector &resultY, Grid &grid, std::vector<int> &curveIndices,
                                          std::vector< std::vector<Intersection> > &mapCurveToConstraints);

    void optimizeImplicitFastWithWeights(Grid &grid, int numberOfVectorFields,
                                         std::vector<PolygonalPath> curves,
                                         std::vector< std::pair<Vector*, Vector*> >& vectorFields,
                                         unsigned short *mapCurveToVectorField,
                                         float *mapCurveToError,
                                         float smoothnessWeight = 0.5);


    //AFKM
public:
    //
    static void AFKMmultiplyByA(const Vector& x, Vector &resultX, Vector &diagATA, AFKMProblemSettings &prob);
    static void optimizeAFKM(Grid &grid, int numberOfAttributeFields, int numberOfAttributes,
                             std::vector<AttributeTrajectory> curves,
                             std::vector< std::vector<Vector > >& attributeFields,
                             unsigned short *mapCurveToVectorField,
                             float *mapCurveToError,
                             double &finalTotalError,
                             int    &totalIterations,
                             std::vector<double> &errorPerAttribute,
                             float smoothnessWeight = 0.5,
                             int maxIterations = 100,
                             AFKMInitializationMethod initializationMethod = FIRST_RANDOM_WORST_MATCH);

    static double AFKMgetTotalError(const std::vector<AFKMCurveDescription> &curves,
                                    const std::vector< std::vector<Vector> > &attributeFields,
                                    const unsigned short *mapCurveToVectorField,
                                    float totalCurveLength,
                                    float smoothnessWeight,
                                    const Grid &grid,
                                    std::vector<double>& totalErrorPerAttribute);
    static double AFKMcomputeErrorImplicit(const std::vector<Vector> &attributeField,
                                           double totalTimeSpan,
                                           const AFKMCurveDescription &curve);
    static double AFKMcomputeErrorImplicit(const std::vector<Vector> &attributeField,
                                           double totalTimeSpan,
                                           const AFKMCurveDescription &curve,
                                           std::vector<double> &errorPerAttrib);
    static void AFKMsetConstraints(std::vector<AFKMCurveDescription> &curve_descriptions,
                                   float &totalCurveLength,
                                   std::vector<AttributeTrajectory> &curves,
                                   const Grid &grid);
    static void AFKMOptimizeAssignments(int &total_change,
                                        double &totalError,
                                        unsigned short *mapCurveToField,
                                        std::vector<std::vector<int> > &mapFieldCurves,
                                        float *mapCurveToError,
                                        const std::vector< std::vector<Vector> > &attributeFields,
                                        const std::vector<AFKMCurveDescription> &curves,
                                        float totalTimeSpan,
                                        float smoothnessWeight,
                                        Grid &grid);

    static void AFKMoptimizeAllAttributeFields(std::vector< std::vector<Vector> > &attributeFields,
                                               Grid &grid,
                                               const std::vector<std::vector<int> > &mapFieldCurves,
                                               const std::vector<AFKMCurveDescription> &curves,
                                               int noOfAttributes,
                                               float totalTimeSpan,
                                               float smoothnessWeight);

    static void AFKMoptimizeField(Grid &grid,
                                  const std::vector<int>& curveIndices,
                                  const std::vector<AFKMCurveDescription> &curve_descriptions,                                  
                                  float totalTimeSpan,
                                  float smoothnessWeight,
                                  int noAttributes,
                                  std::vector<Vector>& initialGuess);

    static void AFKMcg_solve(AFKMProblemSettings &prob,
                             const Vector &b,
                             Vector &x);
    static void AFKMrepopulateEmptyCluster(std::vector<std::vector<int> > &mapFieldCurves,
                                           unsigned short *mapCurveToField,
                                           std::vector< std::vector<Vector> >&);

    static std::pair< std::vector<int>, std::vector< std::vector<int> > > AFKMinitialize(Grid &grid,
                                                                                        int numberOfFields,
                                                                                        int numberOfAttributes,
                                                                                        const std::vector<AFKMCurveDescription> &curves,
                                                                                        float totalTimeSpan,
                                                                                        float smoothnessWeight,
                                                                                        AFKMInitializationMethod initializationMethod);
};

#endif // OPTIMIZER_H
