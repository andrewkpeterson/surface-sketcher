#include "DirectionFieldInitializer.h"
#include "Mesh.h"
#include <Eigen/Sparse>
#include <complex>
#include <map>
#include <cmath>
#include "Sketch.h"

void DirectionFieldInitializer::addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                                         std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b) {

    // here, we are encouraging the a and b variables of the face f to have values such that
    // the direction u is one of the roots of the polynomial that the coefficients a and b define

    // this will encourage the direction u to be one of the four directions in the resulting
    // direction field
    std::complex<double> u(dir.x(), dir.y()); // dir is a unit vector
    std::complex<double> u2 = std::pow(u, 2.0);
    std::complex<double> u4 = std::pow(u, 4.0);
    double cR = u2.real();
    double cI = u2.imag();
    double dR = u4.real();
    double dI = u4.imag();

    {
        double val1 = mult * (cR * cR + cI * cI);
        std::pair<int, int> p1(f->index * 4, f->index * 4);
        addToSparseMap(p1, val1, m);

        double val2 = -mult * cR;
        std::pair<int, int> p2(f->index * 4, (f->index * 4) + 2);
        addToSparseMap(p2, val2, m);

        double val3 = -mult * cI;
        std::pair<int, int> p3(f->index * 4, (f->index * 4) + 3);
        addToSparseMap(p3, val3, m);

        double val4 = -mult * (dR * cR + dI * cI);

        // it is important that we add to the entry of b, rather than just assign to it, otherwise, the
        // constraints will not work correctly
        b(f->index * 4) += val4;
    }

    {
        double val1 = mult * (cR * cR + cI * cI);
        std::pair<int, int> p1((f->index * 4) + 1, (f->index * 4) + 1);
        addToSparseMap(p1, val1, m);

        double val2 = mult * cI;
        std::pair<int, int> p2((f->index * 4) + 1, (f->index * 4) + 2);
        addToSparseMap(p2, val2, m);

        double val3 = -mult * cR;
        std::pair<int, int> p3((f->index * 4) + 1, (f->index * 4) + 3);
        addToSparseMap(p3, val3, m);

        double val4 = mult * (cI * dR - cR * dI);
        b((f->index * 4) + 1) += val4;
    }

    {
        double val1 = -mult * cR;
        std::pair<int, int> p1((f->index * 4) + 2, (f->index * 4));
        addToSparseMap(p1, val1, m);

        double val2 = mult * cI;
        std::pair<int, int> p2((f->index * 4) + 2, (f->index * 4) + 1);
        addToSparseMap(p2, val2, m);

        double val3 = mult;
        std::pair<int, int> p3((f->index * 4) + 2, (f->index * 4) + 2);
        addToSparseMap(p3, val3, m);

        double val4 = mult * dR;
        b((f->index * 4) + 2) += val4;
    }

    {
        double val1 = -mult * cI;
        std::pair<int, int> p1((f->index * 4) + 3, (f->index * 4));
        addToSparseMap(p1, val1, m);

        double val2 = -mult * cR;
        std::pair<int, int> p2((f->index * 4) + 3, (f->index * 4) + 1);
        addToSparseMap(p2, val2, m);

        double val3 = mult;
        std::pair<int, int> p3((f->index * 4) + 3, (f->index * 4) + 3);
        addToSparseMap(p3, val3, m);

        double val4 = mult * dI;
        b((f->index * 4) + 3) += val4;
    }

}

void DirectionFieldInitializer::addCoefficientsForESmooth(Face *f, Face *g, double val, std::map<std::pair<int,int>, double> &m) {

    // here, we are punishing differences in the polyvector representations of neighboring triangles f and g.
    // it makes sense that punishing differences of the polyvector representations of neighboring triangles results
    // in a smoother direction field.
    for (int i = 0; i < 4; i++) {
        std::pair<int, int> p1((f->index * 4) + i, (f->index * 4) + i);
        addToSparseMap(p1, val, m);

        std::pair<int, int> p2((f->index * 4) + i, (g->index * 4) + i);
        addToSparseMap(p2, -val, m);

        std::pair<int, int> p3((g->index * 4) + i, (g->index * 4) + i);
        addToSparseMap(p3, val, m);

        std::pair<int, int> p4((g->index * 4) + i, (f->index * 4) + i);
        addToSparseMap(p4, -val, m);
    }
}

void DirectionFieldInitializer::initializeDirectionFieldFromSolutionVector(Mesh &m, Eigen::VectorXd &x) {
    /*
    m.forEachTriangle([&](std::shared_ptr<Face> f) {
        std::complex<double> a(x(f->index * 4), x((f->index * 4) + 1));
        std::complex<double> b(x((f->index * 4) + 2), x((f->index * 4) + 3));

        std::complex<double> z1 = 1.0 / 2.0 * (a + std::sqrt(a * a - 4.0 * b));
        std::complex<double> z2 = 1.0 / 2.0 * (a - std::sqrt(a * a - 4.0 * b));

        // it should not matter whether the positive or negative square root is taken,
        // because the positive and negative square roots are complex numbers related
        // by a rotation by pi radians, both of which are the directions we want in the
        // 4 direction field
        std::complex<double> u_complex = std::sqrt(z1);
        std::complex<double> v_complex = std::sqrt(z2);
        //std::complex<double> v_complex = std::sqrt(b / z); // we can also solve by v by taking sqrt(a - z), and we will get the same thing

        Eigen::Vector2f u(u_complex.real(), u_complex.imag());
        Eigen::Vector2f v(v_complex.real(), v_complex.imag());

        f->u = u.normalized();
        f->v = v.normalized();
    });
    */

    // go through all of the a values to find the smallest one, indicating the
    // (u,v) pair that is most perpendicular
    int index_with_smallest_a = 0;
    double smallest_a = INFINITY;
    for (int i = 0; i < x.rows(); i += 4) {
        double num = x(i) * x(i) + x(i + 1) * x(i + 1);
        if (num < smallest_a) {
            index_with_smallest_a = i / 4;
            smallest_a = num;
        }
    }

    std::set<Face*> visited_faces;
    initializeDirectionFieldFromSolutionVectorHelper(m, x, visited_faces, index_with_smallest_a, true);
}

void DirectionFieldInitializer::initializeDirectionFieldFromSolutionVectorHelper(Mesh &m, Eigen::VectorXd &x,
                                                                                 std::set<Face*> &visited_faces, int index, bool first) {
    auto f = m.getFace(index);
    visited_faces.insert(f.get());

    std::complex<double> a(x(index * 4), x((index * 4) + 1));
    std::complex<double> b(x((index * 4) + 2), x((index * 4) + 3));
    std::complex<double> z1 = 1.0 / 2.0 * (a + std::sqrt(a * a - 4.0 * b));
    std::complex<double> z2 = 1.0 / 2.0 * (a - std::sqrt(a * a - 4.0 * b));
    std::complex<double> dir1_complex = std::sqrt(z1);
    std::complex<double> dir2_complex = std::sqrt(z2);
    Eigen::Vector2f dir1(dir1_complex.real(), dir1_complex.imag());
    Eigen::Vector2f dir2(dir2_complex.real(), dir2_complex.imag());

    if (first) {
        f->u = dir1;
        f->v = dir2;
    } else {
        Eigen::Vector2f smooth_u = Eigen::Vector2f::Zero();
        int num_initialized_neighbors = 0;
        for (int i = 0; i < f->neighbors.size(); i++) {
            if (f->neighbors[i]->u.norm() > 0) { num_initialized_neighbors++; }
            if (smooth_u.norm() > 0 && smooth_u.dot(f->neighbors[i]->u) < 0) {
                smooth_u += -f->neighbors[i]->u;
            } else {
                smooth_u += f->neighbors[i]->u;
            }
        }
        smooth_u = smooth_u / num_initialized_neighbors;

        Eigen::Vector2f smooth_v = Eigen::Vector2f::Zero();
        for (int i = 0; i < f->neighbors.size(); i++) {
            if (smooth_v.norm() > 0 && smooth_v.dot(f->neighbors[i]->v) < 0) {
                smooth_v += -f->neighbors[i]->v;
            } else {
                smooth_v += f->neighbors[i]->v;
            }
        }
        smooth_v = smooth_v / num_initialized_neighbors;

        Eigen::Vector2f best_u_for_dir1 = (smooth_u.normalized().dot(dir1)) > (-smooth_u.normalized().dot(dir1)) ? smooth_u.normalized() : -smooth_u.normalized();
        Eigen::Vector2f best_v_for_dir1 = (smooth_v.normalized().dot(dir1)) > (-smooth_v.normalized().dot(dir1)) ? smooth_v.normalized() : -smooth_v.normalized();

        Eigen::Vector2f best_u_for_dir2 = (smooth_u.normalized().dot(dir2)) > (-smooth_u.normalized().dot(dir2)) ? smooth_u.normalized() : -smooth_u.normalized();
        Eigen::Vector2f best_v_for_dir2 = (smooth_v.normalized().dot(dir2)) > (-smooth_v.normalized().dot(dir2)) ? smooth_v.normalized() : -smooth_v.normalized();

        bool dir1_assigned_to_v = best_u_for_dir2.dot(dir2) + best_v_for_dir1.dot(dir1) > best_u_for_dir1.dot(dir1) + best_v_for_dir2.dot(dir2);
        //bool dir1_assigned_to_v = smooth_u.dot(dir2) + smooth_v.dot(dir1) > smooth_u.dot(dir1) + smooth_v.dot(dir2);

        if (dir1_assigned_to_v) {
            f->v = dir1;
            f->u = dir2;
        } else {
            f->v = dir2;
            f->u = dir1;
        }
    }

    for (int i = 0; i < f->neighbors.size(); i++) {
        if (visited_faces.find(f->neighbors[i]) == visited_faces.end()) {
            initializeDirectionFieldFromSolutionVectorHelper(m, x, visited_faces, f->neighbors[i]->index, false);
        }
    }
}


void DirectionFieldInitializer::initializeDirectionField(Mesh &mesh, const Sketch &sketch) {

    int num_faces = mesh.getNumTriangles();

    std::vector<Triplet> coefficients;

    SparseMat A(num_faces*4, num_faces*4); // multiply faces by 4 to account for real and imaginary parts of a and b values for each face
    Eigen::VectorXd b = Eigen::VectorXd::Zero(num_faces*4);
    std::map<std::pair<int,int>, double> sparse_map;
    // add coefficients for E_smooth
    mesh.forEachPairOfNeighboringTriangles([&](Face *f, Face *g) { // why are the non-smooth sections always appearing in the same places?
        double efg = Mesh::calcEFGArea(f, g);
        double val = (1.0 / mesh.getTotalArea()) * 2 * efg * OMEGA_S;
        addCoefficientsForESmooth(f, g, val, sparse_map);
    });

    // add coefficients for E_constraint
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        double mult = (1.0 / mesh.getTotalArea()) * 2 * OMEGA_C * f->area;

        if (sketch.checkStrokePoints(f)) {
            for (int i = 0; i < sketch.getStrokePoints(f).size(); i++) {
                Eigen::Vector2f dir = sketch.getStrokePoints(f)[i]->tangent_dir;
                addCoefficientsForEConstraint(f, mult, dir, sparse_map, b);
            }
        }
        if (sketch.checkStrokeLines(f)) {
            for (int i = 0; i < sketch.getStrokeLines(f).size(); i++) {
                Eigen::Vector2f dir = sketch.getStrokeLines(f)[i].dir;
                addCoefficientsForEConstraint(f, mult, dir, sparse_map, b);
            }
        }
    });

    // add coefficients for E_ortho
    // punishing the magnitude of a in the polyvector representation will encourage u and v to be
    // orthogonal. OMEGA_O is small, so using this energy only lightly pushes towards orthogonality.
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        double val = (1.0 / mesh.getTotalArea()) * 2 * OMEGA_O * f->area;
        std::pair<int, int> p1(f->index * 4, f->index * 4);
        addToSparseMap(p1, val, sparse_map);
        std::pair<int, int> p2((f->index * 4) + 1, (f->index * 4) + 1);
        addToSparseMap(p2, val, sparse_map);
    });

    // turn the map into triples and put the triples into the sparse matrix
    for (auto it = sparse_map.begin(); it != sparse_map.end(); it++) {
        coefficients.push_back(Triplet(it->first.first, it->first.second, it->second));
    }
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // solve the system of linear equations
    A.makeCompressed();

    Eigen::SparseLU<SparseMat> solver(A);
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd x = solver.solve(-b); // this is -b rather than b because the solver solves Ax = b, but b was set up to solve Ax + b = 0
    //std::cout << solver.info() << std::endl;
    //std::cout << solver.lastErrorMessage() << std::endl;

    // use the solution to initialize u and v for each face
    initializeDirectionFieldFromSolutionVector(mesh, x);
    checkSmoothness(mesh, x);
}

void DirectionFieldInitializer::checkSmoothness(Mesh &mesh, Eigen::VectorXd &x) {
    double a_max = 0;
    double b_max = 0;
    int num_pairs = 0;
    mesh.forEachPairOfNeighboringTriangles([&](Face *f, Face *g) {
        num_pairs = 0;
        double efg = Mesh::calcEFGArea(f, g);
        double a_new = efg * std::pow(x(f->index * 4) - x(g->index * 4), 2) + std::pow(x(f->index * 4 + 1) - x(g->index * 4 + 1), 2);
        a_max = std::max(a_max, a_new);
        double b_new = efg * std::pow(x(f->index * 4 + 2) - x(g->index * 4 + 2), 2) + std::pow(x(f->index * 4 + 3) - x(g->index * 4 + 3), 2);
        b_max = std::max(b_max, b_new);
        if (b_new > .05 || a_new > .05) {
            std::cout << "*******************" << std::endl;
            std::cout << "a: " << a_new << std::endl;
            std::cout << "b: " << b_new << std::endl;
            std::cout << "x: " << f->centroid.x() << "y: " << f->centroid.y() << std::endl;
            f->discontinuity = true;
            g->discontinuity = true;
        }
    });
    std::cout << "a_max: " << a_max << std::endl;
    std::cout << "b_max: " << b_max << std::endl;
}

void DirectionFieldInitializer::printSparseMatrix(const SparseMat mat) {
    for (int row = 0; row < mat.rows(); row++) {
        for (int col = 0; col < mat.cols(); col++) {
            std::cout << mat.coeff(row, col) << " ";
        }
        std::cout << std::endl;
    }
}
