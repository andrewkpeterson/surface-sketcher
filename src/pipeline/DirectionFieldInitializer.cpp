#include "DirectionFieldInitializer.h"
#include "Mesh.h"
#include <Eigen/Sparse>
#include <complex>
#include <map>
#include <cmath>
#include <queue>
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

    std::map<int, std::pair<Eigen::Vector2f, Eigen::Vector2f>> vectors;

    m.forEachTriangle([&](std::shared_ptr<Face> f) {
        std::complex<double> a(x(f->index * 4), x((f->index * 4) + 1));
        std::complex<double> b(x((f->index * 4) + 2), x((f->index * 4) + 3));

        std::complex<double> z1 = 1.0 / 2.0 * (a + std::sqrt(a * a - 4.0 * b));
        std::complex<double> z2 = 1.0 / 2.0 * (a - std::sqrt(a * a - 4.0 * b));

        // it should not matter whether the positive or negative square root is taken,
        // because the positive and negative square roots are complex numbers related
        // by a rotation by pi radians, both of which are the directions we want in the
        // 4 direction field
        std::complex<double> dir1_complex = std::sqrt(z1);
        std::complex<double> dir2_complex = std::sqrt(z2); // we can also solve for dir2_complex by taking sqrt(a - dir1_complex), and we will get the same thing

        Eigen::Vector2f dir1(dir1_complex.real(), dir1_complex.imag());
        Eigen::Vector2f dir2(dir2_complex.real(), dir2_complex.imag());

        assert(vectors.find(f->index) == vectors.end());
        vectors[f->index] = std::pair(dir1.normalized(), dir2.normalized());
    });

    assignVectors(m, vectors);
    fixSignOfDirectionVectors(m, true);
    fixSignOfDirectionVectors(m, false);
}

void DirectionFieldInitializer::assignVectors(Mesh &m, std::map<int, std::pair<Eigen::Vector2f, Eigen::Vector2f>> vectors) {

    // determine how u, v should be matched to dir1, dir2
    m.forEachTriangle([&](std::shared_ptr<Face> f) {
        int flipped = std::rand() % 2;
        if (flipped) {
            f->u = vectors[f->index].first;
            f->v = vectors[f->index].second;
        } else {
            f->u = vectors[f->index].second;
            f->v = vectors[f->index].first;
        }
    });

    float total_cost = 0;
    m.forEachPairOfNeighboringTriangles([&](Face *f, Face *g) {
        auto u = f->u.dot(g->u) > 0 ? f->u : -f->u;
        auto v = f->v.dot(g->v) > 0 ? f->v : -f->v;
        float u_diff = (u - g->u).norm();
        float v_diff = (v - g->v).norm();
        total_cost += u_diff + v_diff;
    });

    for (int iter = 0; iter < ANNEALING_ITERATIONS; iter++) {
        int idx = std::rand() % m.getNumTriangles();
        Face *f = m.getFace(idx).get();
        float old_cost = 0;
        float new_cost = 0;
        for (int i = 0; i < f->neighbors.size(); i++) {
            auto neighbor = f->neighbors[i];
            auto u_matched_with_neighbor_u = f->u.dot(neighbor->u) > 0 ? f->u : -f->u;
            auto u_matched_with_neighbor_v = f->u.dot(neighbor->v) > 0 ? f->u : -f->u;
            auto v_matched_with_neigbor_v = f->v.dot(neighbor->v) > 0 ? f->v : -f->v;
            auto v_matched_with_neigbor_u = f->v.dot(neighbor->u) > 0 ? f->v : -f->v;
            old_cost += (u_matched_with_neighbor_u - neighbor->u).norm() + (v_matched_with_neigbor_v - neighbor->v).norm();
            new_cost += (v_matched_with_neigbor_u - neighbor->u).norm() + (u_matched_with_neighbor_v - neighbor->v).norm();
        }

        float new_total_cost = total_cost - old_cost + new_cost;

        float temp = INITIAL_TEMPERATURE * std::pow((1 - float(iter) / float(ANNEALING_ITERATIONS)), 2.0);
        float prob = std::exp(-(new_total_cost - total_cost) / temp);
        float rand_num = (float(std::rand()) / float(INT_MAX));
        if (prob >= rand_num) {
            auto saved_v = f->v;
            f->v = f->u;
            f->u = saved_v;
            total_cost = new_total_cost;
        }
    }
}

void DirectionFieldInitializer::fixSignOfDirectionVectors(Mesh &m, bool forV) {
    m.forEachTriangle([&](std::shared_ptr<Face> f) {
        int flipped = std::rand() % 2;
        if (flipped) {
            if (forV) { f->v = -f->v; }
            else { f->u = -f->u; }
        }
    });

    float total_cost = 0;
    m.forEachPairOfNeighboringTriangles([&](Face *f, Face *g) {
        float cost = forV ? f->v.dot(g->v) < 0 : f->u.dot(g->u) < 0;
        total_cost += cost;
    });

    for (int iter = 0; iter < ANNEALING_ITERATIONS_FOR_SIGN; iter++) {
        int idx = std::rand() % m.getNumTriangles();
        Face *f = m.getFace(idx).get();
        float old_cost = 0;
        float new_cost = 0;
        for (int i = 0; i < f->neighbors.size(); i++) {
            auto neighbor = f->neighbors[i];
            old_cost += forV ? f->v.dot(neighbor->v) < 0 : f->u.dot(neighbor->u) < 0;
            new_cost += forV ? -f->v.dot(neighbor->v) < 0 : -f->u.dot(neighbor->u) < 0;
        }

        float new_total_cost = total_cost - old_cost + new_cost;

        float temp = INITIAL_TEMPERATURE * std::pow((1 - float(iter) / float(ANNEALING_ITERATIONS_FOR_SIGN)), 2.0);
        float prob = std::exp(-(new_total_cost - total_cost) / temp);
        float rand_num = (float(std::rand()) / float(INT_MAX));
        if (prob >= rand_num) {
            if (forV) { f->v = -f->v; }
            else { f->u = -f->u; }
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
    mesh.forEachPairOfNeighboringTriangles([&](Face *f, Face *g) {
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
}
