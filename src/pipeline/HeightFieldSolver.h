#ifndef HEIGHTFIELDSOLVER_H
#define HEIGHTFIELDSOLVER_H

#include "Sketch.h"
#include "Mesh.h"

class HeightFieldSolver
{
public:
    static void solveForHeightField(Mesh &mesh, Sketch &sketch);

    static constexpr float PRINCIPLE_CURVATURE_MU = 0;
private:
    static void initializeCurvatureValues(Sketch &sketch);
    static void estimateCurvatureValues(Mesh &mesh, Sketch &sketch);
    static void estimateCurvatureValuesHelper(Mesh &mesh, Sketch &sketch,
                                              std::vector<Sketch::Stroke> &strokes, bool convex);
    static float fitCircle(Mesh &mesh, Sketch &sketch, std::vector<Eigen::Vector2f> points);
    static void minimizeELambda(Mesh &mesh, const Sketch &sketch);
    static void minimizeEMatch(Mesh &mesh, const Sketch &sketch);

    static constexpr int NUM_ITERATIONS = 6;
    static constexpr float INITIAL_CURVATURE_MAGNITUDE = 1;
};

#endif // HEIGHTFIELDSOLVER_H
