#ifndef DIRECTIONFIELDSOLVER_H
#define DIRECTIONFIELDSOLVER_H

#include <map>
#include <complex>
#include <Eigen/Sparse>

class Mesh;
class Sketch;
struct Face;

class DirectionFieldSolver
{
public:
    static void initializeDirectionField(Mesh &mesh, const Sketch &sketch);
    static void optimizeBendFieldEnergy(Mesh &mesh, const Sketch &sketch);

private:
    static void addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                              std::map<std::pair<int,int>, std::complex<double>> &m, Eigen::VectorXcd &b);
    static void addCoefficientsForESmooth(Face *f, Face *g, double mult, std::map<std::pair<int,int>, std::complex<double>> &m);
    static void initializeDirectionFieldFromABVector(Mesh &m, Eigen::VectorXcd &x);

    static constexpr const float OMEGA_C = 10e3;
    static constexpr const float OMEGA_O = 10e-5;
    static constexpr const float SCALE_FACTOR = 100; // factor for scaling A and b so that matrix is not mistaken as non-singular
};

#endif // DIRECTIONFIELDSOLVER_H
