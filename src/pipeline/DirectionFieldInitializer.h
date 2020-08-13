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
    static void printSparseMatrix(const SparseMat);

private:
    static void addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                              std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForESmooth(Face *f, Face *g, double val, std::map<std::pair<int,int>, double> &m);
    static void initializeDirectionFieldFromSolutionVector(Mesh &m, Eigen::VectorXd &x);
    static void initializeDirectionFieldFromSolutionVectorHelper(Mesh &m, Eigen::VectorXd &x, std::set<Face*> &visited_faces, int index, bool first);
    static void checkSmoothness(Mesh &mesh, Eigen::VectorXd &x);
    static inline void addToSparseMap(const std::pair<int, int> &p, double c, std::map<std::pair<int,int>, double> &m) {
        if (m.find(p) == m.end()) {
            m[p] = c;
        } else {
            m[p] += c;
        }
    }

    static constexpr const double OMEGA_S = 1e3; // the smoothness weight was 1 in the paper, but increasing it gives better results
    static constexpr const double OMEGA_C = 1e3; // 1e3
    static constexpr const double OMEGA_O = 1e-5; //1e-5

    static constexpr const double AVERAGE_MULTIPLER = 1.5;
};

#endif // DIRECTIONFIELDSOLVER_H
