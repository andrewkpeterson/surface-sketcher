#include "DirectionFieldOptimizer.h"
#include "Mesh.h"
#include "Sketch.h"

void DirectionFieldOptimizer::optimizeBendFieldEnergy(Mesh &mesh, const Sketch &sketch) {
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        runIterationOfBendFieldEnergyOptimization(mesh, sketch);
    }
}

void DirectionFieldOptimizer::runIterationOfBendFieldEnergyOptimization(Mesh &mesh, const Sketch &sketch) {

    int num_faces = mesh.getNumTriangles();

    // get A and b to calculate new u vectors
    std::vector<Triplet> Au_coefficients;
    SparseMat Au(num_faces*2, num_faces*2);
    Eigen::VectorXd bu = Eigen::VectorXd::Zero(num_faces*2);
    std::map<std::pair<int,int>, double> Au_map;

    addCoefficientsForBendFieldEnergy(mesh, sketch, false, Au_map, bu);
    addCoefficientsForStrokeConstraints(mesh, sketch, false, Au_map, bu);

    for (auto it = Au_map.begin(); it != Au_map.end(); it++) {
        Au_coefficients.push_back(Triplet(it->first.first, it->first.second, it->second));
    }
    Au.setFromTriplets(Au_coefficients.begin(), Au_coefficients.end());

    // solve for u vectors
    Au.makeCompressed();
    Eigen::SparseLU<SparseMat> Usolver(Au);
    //std::cout << Usolver.lastErrorMessage() << std::endl;
    Usolver.analyzePattern(Au);
    //std::cout << Usolver.lastErrorMessage() << std::endl;
    Usolver.factorize(Au);
    //std::cout << "Determinant of matrix: " << Usolver.determinant() << std::endl;
    //std::cout << "Log of determinant of of matrix: " << Usolver.logAbsDeterminant() << std::endl;
    //std::cout << Usolver.lastErrorMessage() << std::endl;
    Eigen::VectorXd xu = Usolver.solve(-bu);
    //std::cout << Usolver.info() << std::endl;
    //std::cout << Usolver.lastErrorMessage() << std::endl;

    // get A and b to calculate new v vectors
    std::vector<Triplet> Av_coefficients;
    SparseMat Av(num_faces*2, num_faces*2);
    Eigen::VectorXd bv = Eigen::VectorXd::Zero(num_faces*2);
    std::map<std::pair<int,int>, double> Av_map;

    addCoefficientsForBendFieldEnergy(mesh, sketch, true, Av_map, bv);
    addCoefficientsForStrokeConstraints(mesh, sketch, true, Av_map, bv);

    for (auto it = Av_map.begin(); it != Av_map.end(); it++) {
        Av_coefficients.push_back(Triplet(it->first.first, it->first.second, it->second));
    }
    Av.setFromTriplets(Av_coefficients.begin(), Av_coefficients.end());

    // solve for v vectors
    Av.makeCompressed();
    Eigen::SparseLU<SparseMat> Vsolver(Av);
    //std::cout << Vsolver.lastErrorMessage() << std::endl;
    Vsolver.analyzePattern(Av);
    //std::cout << Vsolver.lastErrorMessage() << std::endl;
    Vsolver.factorize(Av);
    //std::cout << "Determinant of matrix: " << Vsolver.determinant() << std::endl;
    //std::cout << "Log of determinant of of matrix: " << Vsolver.logAbsDeterminant() << std::endl;
    //std::cout << Vsolver.lastErrorMessage() << std::endl;
    Eigen::VectorXd xv = Vsolver.solve(-bv);
    //std::cout << Vsolver.info() << std::endl;
    //std::cout << Vsolver.lastErrorMessage() << std::endl;

    // normalize u and v vectors and assign new values to mesh
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        Eigen::Vector2f u(xu(2*f->index), xu(2*f->index + 1));
        Eigen::Vector2f v(xv(2*f->index), xv(2*f->index + 1));

        f->u = u.normalized();
        f->v = v.normalized();
    });
}

void DirectionFieldOptimizer::addCoefficientsForStrokeConstraints(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                                  std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b) {
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        if (sketch.checkStrokePoints(f)) {
            for (int i = 0; i < sketch.getStrokePoints(f).size(); i++) {
                Eigen::Vector2f dir = sketch.getStrokePoints(f)[i]->tangent_dir.normalized();
                addCoefficientsForStrokeConstraintsHelper(mesh, sketch, coefficientsForV, m, b, f, dir);
            }
        }

        if (sketch.checkStrokeLines(f)) {
            for (int i = 0; i < sketch.getStrokeLines(f).size(); i++) {
                Eigen::Vector2f dir = sketch.getStrokeLines(f)[i].dir.normalized();
                addCoefficientsForStrokeConstraintsHelper(mesh, sketch, coefficientsForV, m, b, f, dir);
            }
        }
    });
}

void DirectionFieldOptimizer::addCoefficientsForStrokeConstraintsHelper(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                                        std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b,
                                                                        std::shared_ptr<Face> f, const Eigen::Vector2f d) {
    // figure out if this constraint is supposed to be for u or for v
    Eigen::Vector2f best_u = (f->u.normalized().dot(d)) > (-f->u.normalized().dot(d)) ? f->u.normalized() : -f->u.normalized();
    Eigen::Vector2f best_v = (f->v.normalized().dot(d)) > (-f->v.normalized().dot(d)) ? f->v.normalized() : -f->v.normalized();
    //std::cout << (coefficientsForV && (best_v - d).norm() < (best_u - d).norm()) << " " << (!coefficientsForV && (best_u - d).norm() < (best_v - d).norm()) << std::endl;
    bool use_constraint = (coefficientsForV && (best_v.dot(d) > best_u.dot(d))) ||
                          (!coefficientsForV && (best_u.dot(d) > best_v.dot(d)));

    if (use_constraint) {
        int sign = (coefficientsForV && best_v == -f->v) || (!coefficientsForV && best_u == -f->u) ? 1 : -1;
        //std::cout << sign << std::endl;
        float mult = f->area / mesh.getTotalArea();
        std::pair<int, int> x_idx(2*f->index, 2*f->index);
        addToSparseMap(x_idx, mult * STROKE_CONSTRAINT_WEIGHT, m);
        b(2*f->index) = sign * mult * STROKE_CONSTRAINT_WEIGHT * d.x(); // need to take negative of this?

        std::pair<int, int> y_idx(2*f->index + 1, 2*f->index + 1);
        addToSparseMap(y_idx, mult * STROKE_CONSTRAINT_WEIGHT, m);
        b(2*f->index + 1) = sign * mult * STROKE_CONSTRAINT_WEIGHT * d.y(); // need to take negative of this?
    }
}

void DirectionFieldOptimizer::addCoefficientsForBendFieldEnergy(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                                std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b) {
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        Eigen::Vector2f p = f->circumcenter;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 3);
        for (int i = 0; i < f->neighbors.size(); i++) {
            A(0,i) = p.x() - f->neighbors[i]->circumcenter.x();
            A(1,i) = p.y() - f->neighbors[i]->circumcenter.y();
        }

        Eigen::Vector2f w = coefficientsForV ? f->u : f->v; // ******************************* Is this right?
        Eigen::MatrixXd A_pinv = A.completeOrthogonalDecomposition().pseudoInverse();

        float k = w.x() * (A_pinv(0,0) + A_pinv(1,0) + A_pinv(2,0)) + w.y() * (A_pinv(0,1) + A_pinv(1,1) + A_pinv(2,1));
        float k1 = -(w.x() * A_pinv(0,0) + w.y() * A_pinv(0,1));
        float k2 = -(w.x() * A_pinv(1,0) + w.y() * A_pinv(1,1));
        float k3 = -(w.x() * A_pinv(2,0) + w.y() * A_pinv(2,1));

        //std::cout << "x: "<<  w.x() << " " << "y: "<<  w.y() << std::endl;

        assert((k1 == 0 && f->neighbors.size() < 1) || f->neighbors.size() >= 1);
        assert((k2 == 0 && f->neighbors.size() < 2) || f->neighbors.size() >= 2);
        assert((k3 == 0 && f->neighbors.size() < 3) || f->neighbors.size() >= 3);

        float mult = f->area / mesh.getTotalArea();

        {
            int idx = 2*f->index;
            int idx1 = 2*(f->neighbors.size() > 0 ? f->neighbors[0]->index : 0);
            int idx2 = 2*(f->neighbors.size() > 1 ? f->neighbors[1]->index : 0);
            int idx3 = 2*(f->neighbors.size() > 2 ? f->neighbors[2]->index : 0);

            {
                std::pair<int, int> p1(idx, idx);
                addToSparseMap(p1, k * k * mult, m);
                std::pair<int, int> p2(idx, idx1);
                addToSparseMap(p2, k * k1 * mult, m);
                std::pair<int, int> p3(idx, idx2);
                addToSparseMap(p3, k * k2 * mult, m);
                std::pair<int, int> p4(idx, idx3);
                addToSparseMap(p4, k * k3 * mult, m);
            }

            {
                std::pair<int, int> p1(idx1, idx);
                addToSparseMap(p1, k1 * k * mult, m);
                std::pair<int, int> p2(idx1, idx1);
                addToSparseMap(p2, k1 * k1 * mult, m);
                std::pair<int, int> p3(idx1, idx2);
                addToSparseMap(p3, k1 * k2 * mult, m);
                std::pair<int, int> p4(idx1, idx3);
                addToSparseMap(p4, k1 * k3 * mult, m);
            }

            {
                std::pair<int, int> p1(idx2, idx);
                addToSparseMap(p1, k2 * k * mult, m);
                std::pair<int, int> p2(idx2, idx1);
                addToSparseMap(p2, k2 * k1 * mult, m);
                std::pair<int, int> p3(idx2, idx2);
                addToSparseMap(p3, k2 * k2 * mult, m);
                std::pair<int, int> p4(idx2, idx3);
                addToSparseMap(p4, k2 * k3 * mult, m);
            }

            {
                std::pair<int, int> p1(idx3, idx);
                addToSparseMap(p1, k3 * k * mult, m);
                std::pair<int, int> p2(idx3, idx1);
                addToSparseMap(p2, k3 * k1 * mult, m);
                std::pair<int, int> p3(idx3, idx2);
                addToSparseMap(p3, k3 * k2 * mult, m);
                std::pair<int, int> p4(idx3, idx3);
                addToSparseMap(p4, k3 * k3 * mult, m);
            }
        }

        {
            int idx = 2*f->index + 1;
            int idx1 = 2*(f->neighbors.size() > 0 ? f->neighbors[0]->index : 0) + 1;
            int idx2 = 2*(f->neighbors.size() > 1 ? f->neighbors[1]->index : 0) + 1;
            int idx3 = 2*(f->neighbors.size() > 2 ? f->neighbors[2]->index : 0) + 1;

            {
                std::pair<int, int> p1(idx, idx);
                addToSparseMap(p1, k * k * mult, m);
                std::pair<int, int> p2(idx, idx1);
                addToSparseMap(p2, k * k1 * mult, m);
                std::pair<int, int> p3(idx, idx2);
                addToSparseMap(p3, k * k2 * mult, m);
                std::pair<int, int> p4(idx, idx3);
                addToSparseMap(p4, k * k3 * mult, m);
            }

            {
                std::pair<int, int> p1(idx1, idx);
                addToSparseMap(p1, k1 * k * mult, m);
                std::pair<int, int> p2(idx1, idx1);
                addToSparseMap(p2, k1 * k1 * mult, m);
                std::pair<int, int> p3(idx1, idx2);
                addToSparseMap(p3, k1 * k2 * mult, m);
                std::pair<int, int> p4(idx1, idx3);
                addToSparseMap(p4, k1 * k3 * mult, m);
            }

            {
                std::pair<int, int> p1(idx2, idx);
                addToSparseMap(p1, k2 * k * mult, m);
                std::pair<int, int> p2(idx2, idx1);
                addToSparseMap(p2, k2 * k1 * mult, m);
                std::pair<int, int> p3(idx2, idx2);
                addToSparseMap(p3, k2 * k2 * mult, m);
                std::pair<int, int> p4(idx2, idx3);
                addToSparseMap(p4, k2 * k3 * mult, m);
            }

            {
                std::pair<int, int> p1(idx3, idx);
                addToSparseMap(p1, k3 * k * mult, m);
                std::pair<int, int> p2(idx3, idx1);
                addToSparseMap(p2, k3 * k1 * mult, m);
                std::pair<int, int> p3(idx3, idx2);
                addToSparseMap(p3, k3 * k2 * mult, m);
                std::pair<int, int> p4(idx3, idx3);
                addToSparseMap(p4, k3 * k3 * mult, m);
            }
        }

    });
}
