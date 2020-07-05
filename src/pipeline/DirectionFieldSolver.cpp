#include "DirectionFieldSolver.h"
#include "Mesh.h"
#include <Eigen/Sparse>
#include <complex>
#include <map>
#include <cmath>
#include "Sketch.h"

inline void addToSparseMap(const std::pair<int, int> &p, std::complex<double> c, std::map<std::pair<int,int>, std::complex<double>> &m) {
    if (m.find(p) == m.end()) {
        m[p] = c;
    } else {
        m[p] += c;
    }
}

void DirectionFieldSolver::addCoefficientsForEConstraint(std::shared_ptr<Face> f, double mult, Eigen::Vector2f dir,
                                                         std::map<std::pair<int,int>, std::complex<double>> &m, Eigen::VectorXcd &b) {
    std::complex<double> u(dir.x(), dir.y());

    std::complex<double> val1 = mult * std::pow(u, 4.0);
    std::pair<int, int> p1(f->index * 2, f->index * 2);
    addToSparseMap(p1, val1, m);

    std::complex<double> val2 = -mult * std::pow(u, 2.0f);
    std::pair<int, int> p2(f->index * 2, (f->index * 2) + 1);
    addToSparseMap(p2, val2, m);

    std::complex<double> val3 = -mult * std::pow(u, 6.0f);
    b(f->index * 2) = val3;

    std::complex<double> val4 = std::complex<float>(mult, 0);
    std::pair<int, int> p4((f->index * 2) + 1, (f->index) * 2 + 1);
    addToSparseMap(p4, val4, m);

    std::complex<double> val5 = -mult * std::pow(u, 2.0f);
    std::pair<int, int> p5((f->index * 2) + 1, (f->index * 2));
    addToSparseMap(p5, val5, m);

    std::complex<double> val6 = mult * std::pow(u, 4.0f);
    b((f->index * 2) + 1) = val6;
}

void DirectionFieldSolver::addCoefficientsForESmooth(Face *f, Face *g, double mult, std::map<std::pair<int,int>, std::complex<double>> &m) {
    std::complex<double> c(mult, 0);

    std::pair<int, int> p1(f->index * 2, f->index * 2);
    addToSparseMap(p1, c, m);

    std::pair<int, int> p2(f->index * 2, g->index * 2);
    addToSparseMap(p2, -c, m);

    std::pair<int, int> p3((f->index * 2) + 1, (f->index * 2) + 1);
    addToSparseMap(p3, c, m);

    std::pair<int, int> p4((f->index * 2) + 1, (g->index * 2) + 1);
    addToSparseMap(p4, -c, m);

    std::pair<int, int> p5(g->index * 2, g->index * 2);
    addToSparseMap(p5, c, m);

    std::pair<int, int> p6(g->index * 2, f->index * 2);
    addToSparseMap(p6, -c, m);

    std::pair<int, int> p7((g->index * 2) + 1, (g->index * 2) + 1);
    addToSparseMap(p7, c, m);

    std::pair<int, int> p8((g->index * 2) + 1, (f->index * 2) + 1);
    addToSparseMap(p8, -c, m);
}

void DirectionFieldSolver::initializeDirectionFieldFromABVector(Mesh &m, Eigen::VectorXcd &x) {
    m.forEachTriangle([&](std::shared_ptr<Face> f) {
        std::complex<double> a = x(f->index * 2);
        std::complex<double> b = x((f->index * 2) + 1);

        std::complex<double> z = 1.0 / 2.0 * (a + std::sqrt(a * a - 4.0 * b));
        std::complex<double> u_complex = std::sqrt(z); // does it matter if the positive or negative square root is taken?
        std::complex<double> v_complex = std::sqrt(b / z);

        Eigen::Vector2f u(u_complex.real(), u_complex.imag());
        Eigen::Vector2f v(v_complex.real(), v_complex.imag());

        f->u = u.normalized();
        f->v = v.normalized();
    });
}

void DirectionFieldSolver::initializeDirectionField(Mesh &mesh, const Sketch &sketch) {
    using SparseMat = Eigen::SparseMatrix<std::complex<double>>;
    using CTriplet = Eigen::Triplet<std::complex<double>>;

    int num_faces = mesh.getNumTriangles();

    std::vector<CTriplet> coefficients;

    SparseMat A(num_faces*2, num_faces*2); // multiply faces by 2 to account for a and b values for each face
    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(num_faces*2);
    std::map<std::pair<int,int>, std::complex<double>> sparse_map;

    // add coefficients for E_smooth
    mesh.forEachPairOfNeighboringTriangles([&](Face *f, Face *g) {
        double efg = Mesh::calcEFGArea(f, g);
        double mult = (1.0 / mesh.getTotalArea()) * 2 * efg * SCALE_FACTOR;
        addCoefficientsForESmooth(f, g, mult, sparse_map);
    });

    // add coefficients for E_constraint
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        double mult = (1.0 / mesh.getTotalArea()) * 2 * OMEGA_C * f->area * SCALE_FACTOR;

        if (sketch.checkStrokePoints(f)) {
            for (int i = 0; i < sketch.getConstStrokePoints(f).size(); i++) {
                Eigen::Vector2f dir = sketch.getConstStrokePoints(f)[i]->tangent_dir;
                addCoefficientsForEConstraint(f, mult, dir, sparse_map, b);
            }
        }

        if (sketch.checkStrokeLines(f)) {
            for (int i = 0; i < sketch.getConstStrokeLines(f).size(); i++) {
                Eigen::Vector2f dir = sketch.getConstStrokeLines(f)[i].dir;
                addCoefficientsForEConstraint(f, mult, dir, sparse_map, b);
            }
        }
    });

    // add coefficients for E_ortho
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        int idx = f->index * 2;
        double val = (1.0 / mesh.getTotalArea()) * 2 * OMEGA_O * f->area * SCALE_FACTOR;
        std::complex<double> c(val,0);
        std::pair<int, int> p(idx, idx);
        addToSparseMap(p, c, sparse_map);
    });

    // turn the map into triples and put the triples into the sparse matrix
    for (auto it = sparse_map.begin(); it != sparse_map.end(); it++) {
        coefficients.push_back(CTriplet(it->first.first, it->first.second, it->second));
    }
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // solve the system of linear equations
    A.makeCompressed();

    //Eigen::LeastSquaresConjugateGradient<SparseMat> solver(A);
    //Eigen::BiCGSTAB<SparseMat> solver(A);
    //Eigen::SparseQR<SparseMat, Eigen::COLAMDOrdering<int>> solver(A);

    Eigen::SparseLU<SparseMat> solver(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    solver.analyzePattern(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    solver.factorize(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    Eigen::VectorXcd x = solver.solve(-b);
    std::cout << solver.lastErrorMessage() << std::endl;
    // use the solution to initialize u and v for each face
    initializeDirectionFieldFromABVector(mesh, x);
}

void DirectionFieldSolver::optimizeBendFieldEnergy(Mesh &mesh, const Sketch &sketch) {

}
