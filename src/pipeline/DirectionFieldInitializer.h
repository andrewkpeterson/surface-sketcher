#ifndef DIRECTIONFIELDSOLVER_H
#define DIRECTIONFIELDSOLVER_H

#include <map>
#include <complex>
#include <Eigen/Sparse>

class Mesh;
class Sketch;
struct Face;

using SparseMat = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

class DirectionFieldInitializer
{
public:
    static void initializeDirectionField(Mesh &mesh, const Sketch &sketch);
    static void printSparseMatrix(const SparseMat);

private:
    static void addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                              std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForESmooth(Face *f, Face *g, double val, std::map<std::pair<int,int>, double> &m);
    static void initializeDirectionFieldFromABVector(Mesh &m, Eigen::VectorXd &x);
    static inline void addToSparseMap(const std::pair<int, int> &p, double c, std::map<std::pair<int,int>, double> &m) {
        if (m.find(p) == m.end()) {
            m[p] = c;
        } else {
            m[p] += c;
        }
    }

    static constexpr const double OMEGA_S = 10;
    static constexpr const double OMEGA_C = 1e3; // 1e3
    static constexpr const double OMEGA_O = 1e-5; //1e-5
};

#endif // DIRECTIONFIELDSOLVER_H
