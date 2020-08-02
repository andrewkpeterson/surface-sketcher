#ifndef DIRECTIONFIELDOPTIMIZER_H
#define DIRECTIONFIELDOPTIMIZER_H

#include "DirectionFieldInitializer.h"
#include <iostream>

class DirectionFieldOptimizer
{
public:
    static void optimizeBendFieldEnergy(Mesh &mesh, const Sketch &sketch);

private:
    static void runIterationOfBendFieldEnergyOptimization(Mesh &mesh, const Sketch &sketch);
    static void addCoefficientsForStrokeConstraints(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                    std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForBendFieldEnergy(Mesh &mesh, const Sketch &sketch, bool coefficientsForV, std::map<std::pair<int,int>,
                                                  double> &m, Eigen::VectorXd &b);
    static void addCoefficientsForStrokeConstraintsHelper(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                          std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b,
                                                          std::shared_ptr<Face> f, const Eigen::Vector2f d);

    static inline void addToSparseMap(const std::pair<int, int> &p, double c, std::map<std::pair<int,int>, double> &m) {
        if (m.find(p) == m.end()) {
            m[p] = c;
        } else {
            m[p] += c;
        }
    }

    static void printSparseMatrix(const SparseMat &mat) {
        for (int row = 0; row < mat.rows(); row++) {
            for (int col = 0; col < mat.cols(); col++) {
                std::cout << mat.coeff(row, col) << " ";
            }
            std::cout << std::endl;
        }
    }

    static void printSparseMatrixCol(const SparseMat &mat, int col) {
        for (int row = 0; row < mat.rows(); row++) {
            std::cout << mat.coeff(row, col) << " ";
        }
    }

    static constexpr const float STROKE_CONTRAINT_WEIGHT = .1;
    static constexpr const int NUM_ITERATIONS = 5;
};

#endif // DIRECTIONFIELDOPTIMIZER_H
