#include "DirectionFieldOptimizer.h"
#include "Mesh.h"
#include "Sketch.h"

void DirectionFieldOptimizer::optimizeBendFieldEnergy(Mesh &mesh, Sketch &sketch) {
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
    Usolver.analyzePattern(Au);
    Usolver.factorize(Au);
    Eigen::VectorXd xu = Usolver.solve(-bu); // this is -bu rather than bu because the solver solves Ax = bu, but b was set up to solve Ax + bu = 0
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
    Vsolver.analyzePattern(Av);
    Vsolver.factorize(Av);
    Eigen::VectorXd xv = Vsolver.solve(-bv); // this is -bv rather than bv because the solver solves Ax = bv, but b was set up to solve Ax + bv = 0
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
        std::set<int> strokes;
        if (sketch.checkStrokePoints(f)) {
            for (int i = 0; i < sketch.getStrokePoints(f).size(); i++) {
                if (strokes.find(sketch.getStrokePoints(f)[i]->strokeId) == strokes.end()) {
                    auto point = sketch.getStrokePoints(f)[i];
                    Eigen::Vector2f dir = sketch.getStrokePoints(f)[i]->tangent_dir.normalized();
                    addCoefficientsForStrokeConstraintsHelper(mesh, sketch, coefficientsForV, m, b, f, dir);
                }
                strokes.insert(sketch.getStrokePoints(f)[i]->strokeId);
            }
        }

        if (sketch.checkStrokeLines(f)) {
            for (int i = 0; i < sketch.getStrokeLines(f).size(); i++) {
                if (strokes.find(sketch.getStrokeLines(f)[i].points.first->strokeId) == strokes.end()) {
                    auto point = sketch.getStrokeLines(f)[i];
                    Eigen::Vector2f dir = sketch.getStrokeLines(f)[i].dir.normalized();
                    addCoefficientsForStrokeConstraintsHelper(mesh, sketch, coefficientsForV, m, b, f, dir);
                }
                strokes.insert(sketch.getStrokeLines(f)[i].points.first->strokeId);
            }
        }
    });
}

void DirectionFieldOptimizer::addCoefficientsForStrokeConstraintsHelper(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                                        std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b,
                                                                        std::shared_ptr<Face> f, const Eigen::Vector2f d) {
    // We need to figure out if this constraint is supposed to be for u or for v.
    // The constraint is supposed to be for u if f->u or -f->u is closer to the
    // constraint than both of f->v and -f->v.
    Eigen::Vector2f smooth_u(0,0);
    Eigen::Vector2f smooth_v(0,0);
    for (int i = 0; i < f->neighbors.size(); i++) {
        smooth_u += f->neighbors[i]->u;
        smooth_v += f->neighbors[i]->v;
    }
    smooth_u = smooth_u / (f->neighbors.size() + 1);
    smooth_v = smooth_v / (f->neighbors.size() + 1);
    Eigen::Vector2f best_u = (smooth_u.normalized().dot(d) > 0) ? smooth_u.normalized() : -smooth_u.normalized();
    Eigen::Vector2f best_v = (smooth_v.normalized().dot(d) > 0) ? smooth_v.normalized() : -smooth_v.normalized();
    bool use_constraint = (coefficientsForV && (best_v.dot(d) > best_u.dot(d))) ||
                          (!coefficientsForV && (best_u.dot(d) > best_v.dot(d)));


    if (use_constraint) {
        if (coefficientsForV && std::abs(d.x()) > std::abs(d.y())) {
            std::cout << "stop" << std::endl;
        }
        // Here, we are creating an equation of the form Ax + b = 0. Note that the Eigen solver actually solves Wx = z, so
        // we pass -b into the solver. The modified constraint is the constraint or the negative of the constraint,
        // whichever is closer to the vector being constrained. To create the equation of the form Ax + b, we want to add
        // the negative of the modified constraint to b at the correct indices.
        Eigen::Vector2f modified_constraint = ((coefficientsForV && f->v.dot(d) > 0) || (!coefficientsForV && f->u.dot(d) > 0)) ? d.normalized() : -d.normalized();
        float mult = 2 * f->area / mesh.getTotalArea() * STROKE_CONSTRAINT_WEIGHT;
        std::pair<int, int> x_idx(2*f->index, 2*f->index);
        addToSparseMap(x_idx, mult, m);
        b(2*f->index) += -mult * modified_constraint.x();

        std::pair<int, int> y_idx(2*f->index + 1, 2*f->index + 1);
        addToSparseMap(y_idx, mult, m);
        b(2*f->index + 1) += -mult * modified_constraint.y();
    }
}

void DirectionFieldOptimizer::addCoefficientsForBendFieldEnergy(Mesh &mesh, const Sketch &sketch, bool coefficientsForV,
                                                                std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b) {
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        Eigen::Vector2f p = f->centroid;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 3);
        // TODO: figure out if circumcenter is better!
        for (int i = 0; i < f->neighbors.size(); i++) {
            A(0,i) = f->neighbors[i]->centroid.x() - p.x();
            A(1,i) = f->neighbors[i]->centroid.y() - p.y();
        }

        // If we are optimizing the v direction field, we multiply the jacobian of v with
        // the u of the same face and then take the norm to get the bendfield energy
        Eigen::Vector2f w = coefficientsForV ? f->u : f->v;
        Eigen::MatrixXd A_pinv = A.completeOrthogonalDecomposition().pseudoInverse();

        float k = w.x() * (A_pinv(0,0) + A_pinv(1,0) + A_pinv(2,0)) + w.y() * (A_pinv(0,1) + A_pinv(1,1) + A_pinv(2,1));
        float k1 = -(w.x() * A_pinv(0,0) + w.y() * A_pinv(0,1));
        float k2 = -(w.x() * A_pinv(1,0) + w.y() * A_pinv(1,1));
        float k3 = -(w.x() * A_pinv(2,0) + w.y() * A_pinv(2,1));

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
