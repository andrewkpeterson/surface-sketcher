#ifndef DIRECTIONFIELDSOLVER_H
#define DIRECTIONFIELDSOLVER_H

#include <map>
#include <set>
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

private:
    static void addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                              std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForESmooth(Face *f, Face *g, double val, std::map<std::pair<int,int>, double> &m);

    static void initializeDirectionFieldFromSolutionVector(Mesh &m, Eigen::VectorXd &x);
    static void assignVectors(Mesh &m, std::map<int, std::pair<Eigen::Vector2f, Eigen::Vector2f>>);
    static void fixSignOfDirectionVectors(Mesh &m, bool forV);
    static inline void addToSparseMap(const std::pair<int, int> &p, double c, std::map<std::pair<int,int>, double> &m) {
        if (m.find(p) == m.end()) {
            m[p] = c;
        } else {
            m[p] += c;
        }
    }

    static constexpr const double OMEGA_S = 1e3; // the smoothness weight was 1 in the paper, but increasing it to 1e3 up to 8e3 gives better results
    static constexpr const double OMEGA_C = 1e3; // 1e3
    static constexpr const double OMEGA_O = 1e-5; //1e-5

    static constexpr const int ANNEALING_ITERATIONS = 2e8; // 2e8 for .2 fineness, 2e8 for .1 fineness, 6e8 for .05 fineness
    static constexpr const int ANNEALING_ITERATIONS_FOR_SIGN = 2e8; // 2e8 for .2 fineness, 2e8 for .1 fineness, 6e8 for .05 fineness
    static constexpr const float INITIAL_TEMPERATURE = 7; // 5 for .2 fineness, 5 for .1 fineness, 5 for .05 fineness
};

#endif // DIRECTIONFIELDSOLVER_H
