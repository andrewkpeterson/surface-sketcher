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
        std::complex<double> dir2_complex = std::sqrt(z2);
        //std::complex<double> v_complex = std::sqrt(b / z); // we can also solve by v by taking sqrt(a - z), and we will get the same thing

        Eigen::Vector2f dir1(dir1_complex.real(), dir1_complex.imag());
        Eigen::Vector2f dir2(dir2_complex.real(), dir2_complex.imag());

        assert(vectors.find(f->index) == vectors.end());
        vectors[f->index] = std::pair(dir1.normalized(), dir2.normalized());
    });

    // go through all of the a values to find the smallest one, indicating the
    // (u,v) pair that is most perpendicular
    int most_perpendicular_index = 0;
    double smallest_dot = INFINITY;
    for (auto it = vectors.begin(); it != vectors.end(); it++) {
        double dot = it->second.first.dot(it->second.second);
        if (dot < smallest_dot) {
            most_perpendicular_index = it->first;
            smallest_dot = dot;
        }
    }

    std::set<Face*> visited_faces_initial;
    Eigen::VectorXf initial_assignments = Eigen::VectorXf::Zero(m.getNumTriangles());
    auto f = m.getFace(most_perpendicular_index);
    Eigen::Vector2f dir1 = vectors[most_perpendicular_index].first;
    Eigen::Vector2f dir2 = vectors[most_perpendicular_index].second;
    f->u = dir1;
    f->v = dir2;
    visited_faces_initial.insert(f.get());

    // this is an improvement on depth-first
    initializeDirectionFieldFromSolutionVectorHelper(m, vectors, visited_faces_initial, most_perpendicular_index);

    //make all u vectors flow same way, do same for v vectors
    std::set<Face*> visited_faces_fix_direction;
    fixSignOfDirectionVectors(m, most_perpendicular_index, visited_faces_fix_direction);
}

void DirectionFieldInitializer::fixSignOfDirectionVectors(Mesh &m, int index, std::set<Face*> &visited_faces) {
    auto f = m.getFace(index);
    int negate_u = 0;
    int negate_v = 0;
    int num_neighbors = 0;
    for (int i = 0; i < f->neighbors.size(); i++) {
        if (visited_faces.find(f->neighbors[i]) != visited_faces.end()) {
            num_neighbors++;
            if (f->u.dot(f->neighbors[i]->u) < 0) { negate_u++; }
            if (f->v.dot(f->neighbors[i]->v) < 0) { negate_v++; }
        }
    }
    if (negate_u > num_neighbors) { f->u = -f->u; }
    if (negate_v > num_neighbors) { f->v = -f->v; }
    visited_faces.insert(f.get());
    for (int i = 0; i < f->neighbors.size(); i++) {
        if (visited_faces.find(f->neighbors[i]) == visited_faces.end()) {
            fixSignOfDirectionVectors(m, f->neighbors[i]->index, visited_faces);
        }
    }
}

Eigen::Vector2f DirectionFieldInitializer::averageNearbyVectors(Mesh &m, std::set<Face*> &visited_faces, std::set<Face*> &added_faces, Face *f, int depth_left, bool forV) {
    added_faces.insert(f);
    if (depth_left > 0) {
        Eigen::Vector2f sum(0,0);
        for (int i = 0; i < f->neighbors.size(); i++) {
            if (visited_faces.find(f->neighbors[i]) != visited_faces.end() && added_faces.find(f->neighbors[i]) == added_faces.end()) {
                Eigen::Vector2f new_vec = averageNearbyVectors(m, visited_faces, added_faces, f->neighbors[i], depth_left-1, forV);
                if (sum.dot(new_vec) < 0) {
                    sum += -new_vec;
                } else {
                    sum += new_vec;
                }
            }
        }

        Eigen::Vector2f curr_face_vector = forV ? f->v : f->u;
        curr_face_vector = sum.dot(curr_face_vector) < 0 ? -curr_face_vector : curr_face_vector;

        if (visited_faces.find(f) != visited_faces.end() && (sum + curr_face_vector).norm() > 0 && depth_left < AVERAGING_DEPTH) { return (sum + curr_face_vector); }
        else if (sum.norm() > 0) { return sum; }
        else { return Eigen::Vector2f(0,0); }
    } else {
        Eigen::Vector2f curr_face_vector = forV ? f->v : f->u;
        if (visited_faces.find(f) != visited_faces.end() && curr_face_vector.norm() > 0) { return curr_face_vector; }
        else { return Eigen::Vector2f(0,0); }
    }
}

void DirectionFieldInitializer::initializeDirectionFieldFromSolutionVectorHelper(Mesh &m, std::map<int, std::pair<Eigen::Vector2f, Eigen::Vector2f>> &vectors,
                                                                                        std::set<Face*> &visited_faces, int start_index) {
    std::queue<Face*> queue;
    auto f = m.getFace(start_index);
    visited_faces.insert(f.get());
    for (int i = 0; i < f->neighbors.size(); i++) {
        queue.push(f->neighbors[i]);
    }

    while(!queue.empty()) {
        auto f = queue.front();
        queue.pop();
        visited_faces.insert(f);

        Eigen::Vector2f dir1 = vectors[f->index].first;
        Eigen::Vector2f dir2 = vectors[f->index].second;

        std::set<Face*> added_faces_u;
        Eigen::Vector2f smooth_u = averageNearbyVectors(m, visited_faces, added_faces_u, f, AVERAGING_DEPTH, false).normalized();
        std::set<Face*> added_faces_v;
        Eigen::Vector2f smooth_v = averageNearbyVectors(m, visited_faces, added_faces_v, f, AVERAGING_DEPTH, true).normalized();

        Eigen::Vector2f best_u_for_dir1 = (smooth_u.dot(dir1)) > (-smooth_u.dot(dir1)) ? smooth_u : -smooth_u;
        Eigen::Vector2f best_v_for_dir1 = (smooth_v.dot(dir1)) > (-smooth_v.dot(dir1)) ? smooth_v : -smooth_v;
        Eigen::Vector2f best_u_for_dir2 = (smooth_u.dot(dir2)) > (-smooth_u.dot(dir2)) ? smooth_u : -smooth_u;
        Eigen::Vector2f best_v_for_dir2 = (smooth_v.dot(dir2)) > (-smooth_v.dot(dir2)) ? smooth_v : -smooth_v;

        bool dir1_assigned_to_v = best_u_for_dir2.dot(dir2) + best_v_for_dir1.dot(dir1) > best_u_for_dir1.dot(dir1) + best_v_for_dir2.dot(dir2);

        if (dir1_assigned_to_v) {
            f->v = dir1;
            f->u = dir2;
        } else {
            f->v = dir2;
            f->u = dir1;
        }

        for (int i = 0; i < f->neighbors.size(); i++) {
            if (visited_faces.find(f->neighbors[i]) == visited_faces.end()) {
                queue.push(f->neighbors[i]);
            }
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
}

void DirectionFieldInitializer::printSparseMatrix(const SparseMat mat) {
    for (int row = 0; row < mat.rows(); row++) {
        for (int col = 0; col < mat.cols(); col++) {
            std::cout << mat.coeff(row, col) << " ";
        }
        std::cout << std::endl;
    }
}
