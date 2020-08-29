#ifndef HEIGHTFIELDSOLVER_H
#define HEIGHTFIELDSOLVER_H

#include "Sketch.h"
#include "Mesh.h"
#include "DirectionFieldInitializer.h"

class HeightFieldSolver
{
public:
    static void solveForHeightField(Mesh &mesh, Sketch &sketch);

    // if the radius of the estimated circle is not regualrized, we will underestimate the normal curvature
    // and the surface curvature will vanish
    static constexpr double NORMAL_CURVATURE_MU = 1e-6; //2e-5; // 1e-6 // this was 2e-5 in the paper

private:

    static void initializeCurvatureValues(Sketch &sketch);
    static void estimateCurvatureValues(Mesh &mesh, Sketch &sketch);
    static void estimateCurvatureValuesHelper(Mesh &mesh, Sketch &sketch,
                                              const std::vector<Sketch::Stroke> &strokes, bool convex);
    static void estimateCurvatureValuesHelperMorePrecise(Mesh &mesh, Sketch &sketch,
                                              const std::vector<Sketch::Stroke> &strokes, bool convex);
    static float fitCircle(Mesh &mesh, Sketch &sketch, std::vector<Eigen::Vector2f> points);

    static void minimizeELambda(Mesh &mesh, const Sketch &sketch);
    static void addCoefficientsForELambda(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>,double> &m, bool constraining_lambda_u, bool in_dir_u);
    static void addConstraintsForELambda(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int, int>, double> &m, Eigen::VectorXd &b);
    static void addConstraintsForELambdaHelper(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int, int>, double> &m, std::shared_ptr<Face> f, Eigen::VectorXd &b,
                                               Eigen::Vector2f dir, float curvature_value);

    struct Values {
        float n3;
        float G11;
        float G12;
        float G21;
        float G22;
        float a11;
        float a12;
        float a13;
        float a21;
        float a22;
        float a23;
        float a31;
        float a32;
        float a33;
        std::vector<Face*> faces;
    };

    static void optimizeHeightField(Mesh &mesh, const Sketch &sketch);
    static void optimizeHeightFieldHelper(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void computeAndAddCoefficientsForEMatch(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b,
                                                   std::shared_ptr<Face> f, bool u_direction, bool x_coordinate, HeightFieldSolver::Values &vs);
    static Eigen::Matrix2f calculateGInverse(std::shared_ptr<Face> f);
    static float calcK(Face *f, std::shared_ptr<Vertex> v, int face_idx, Values &vs, bool u_direction, bool x_coordinate); // idx == 0 indicates center
    static float calcdN(Face *f, std::shared_ptr<Vertex> v, int face_idx, Values &vs, int N_idx);
    static float calcdAlpha(Face *f, std::shared_ptr<Vertex> v, int face_idx, Values &vs);
    static float calcdBeta(Face *f, std::shared_ptr<Vertex> v, int face_idx, Values &vs);
    static float calcdGamma(Face *f, std::shared_ptr<Vertex> v, int face_idx, Values &vs);
    static float calcdB(Face *f, std::shared_ptr<Vertex> v, int face_idx, int B_idx);
    static float calcdG(Face *f, std::shared_ptr<Vertex> v, int face_idx, bool deriv_wrt_x);
    static void addCoefficientsToMapAndVectorForEMatch(Face* f, Mesh &mesh, std::vector<float> coefficients, std::vector<std::shared_ptr<Vertex>> vertices,
                                                       std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b, float constraint);
    static void addCoefficientsForEBoundary(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForRegualrityConstraint(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m);
    static void addCoefficientsForRegualrityConstraintHelper(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m,
                                                             std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2);
    static bool faceContainsEdgeVertex(Mesh &mesh, Face *f);
    static std::shared_ptr<Vertex> getEdgeVertexNotInSet(Mesh &mesh, Face *f, const std::set<std::shared_ptr<Vertex> > &set);
    static void addCoefficientsForEContour(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int, int>, double> &m, Eigen::VectorXd &b);

    static inline void addToSparseMap(const std::pair<int, int> &p, double c, std::map<std::pair<int,int>, double> &m) {
        if (m.find(p) == m.end()) {
            m[p] = c;
        } else {
            m[p] += c;
        }
    }

    static constexpr int NUM_ITERATIONS = 20;
    static constexpr float INITIAL_CURVATURE_MAGNITUDE = 1;

    static constexpr float CURVATURE_MAGNITUDE_SMOOTHNESS_BETA = .01;
    static constexpr float CURVATURE_MAGNITUDE_CONSTRAINT_WEIGHT = 1e6;

    static constexpr float BOUNDARY_POSITIONAL_CONSTRAINT_WEIGHT = 10; // omega_0, changes for different cases
    static constexpr float BOUNDARY_REGULARITY_CONSTRAINT_WEIGHT = 1; // omega_1, same for all cases
    static constexpr float BOUNDARY_HEIGHT = 0;

    static constexpr float CONTOUR_CONSTRAINT_WEIGHT = .5;
};

#endif // HEIGHTFIELDSOLVER_H
