#ifndef DIRECTIONFIELDSOLVER_H
#define DIRECTIONFIELDSOLVER_H

#include <map>
#include <complex>
#include <Eigen/Sparse>

class Mesh;
class Sketch;
struct Face;

using SparseMat = Eigen::SparseMatrix<double>;
using CTriplet = Eigen::Triplet<double>;

class DirectionFieldSolver
{
public:
    static void initializeDirectionField(Mesh &mesh, const Sketch &sketch);
    static void optimizeBendFieldEnergy(Mesh &mesh, const Sketch &sketch);
    static void printSparseMatrix(const SparseMat);

private:
    static void addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                              std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForESmooth(Face *f, Face *g, double val, std::map<std::pair<int,int>, double> &m);
    static void initializeDirectionFieldFromABVector(Mesh &m, Eigen::VectorXd &x);

    static constexpr const double OMEGA_C = 1; //10e3;
    static constexpr const double OMEGA_O = 1; //10e-5;
};

#endif // DIRECTIONFIELDSOLVER_H
