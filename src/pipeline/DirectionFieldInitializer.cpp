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
        b(f->index * 4) = val4;
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
        b((f->index * 4) + 1) = val4;
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
        b((f->index * 4) + 2) = val4;
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
        b((f->index * 4) + 3) = val4;
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

void DirectionFieldInitializer::initializeDirectionFieldFromABVector(Mesh &m, Eigen::VectorXd &x) {
    m.forEachTriangle([&](std::shared_ptr<Face> f) {
        std::complex<double> a(x(f->index * 4), x((f->index * 4) + 1));
        std::complex<double> b(x((f->index * 4) + 2), x((f->index * 4) + 3));

        std::complex<double> z = 1.0 / 2.0 * (a + std::sqrt(a * a - 4.0 * b));

        // it should not matter whether the positive or negative square root is taken,
        // because the positive and negative square roots are complex numbers related
        // by a rotation by pi radians, both of which are the directions we want in the
        // 4 direction field
        std::complex<double> u_complex = std::sqrt(z);
        std::complex<double> v_complex = std::sqrt(b / z); // we can also solve by v by taking sqrt(a - z), and we will get the same thing
        std::complex<double> u_complex_normalized = u_complex / (std::sqrt(u_complex.real() * u_complex.real() + u_complex.imag() * u_complex.imag()));
        std::complex<double> v_complex_normalized = v_complex / (std::sqrt(v_complex.real() * v_complex.real() + v_complex.imag() * v_complex.imag()));
        std::cout << (std::pow(u_complex_normalized, 4.0) - std::pow(u_complex_normalized, 2.0) * a + b) << std::endl;
        std::cout << (std::pow(v_complex_normalized, 4.0) - std::pow(v_complex_normalized, 2.0) * a + b) << std::endl;

        Eigen::Vector2f u(u_complex.real(), u_complex.imag());
        Eigen::Vector2f v(v_complex.real(), v_complex.imag());

        f->u = u.normalized();
        f->v = v.normalized();
    });
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

    //printSparseMatrix(A);

    /*
    Eigen::ConjugateGradient<SparseMat> solverCG;
    auto A_transpose = SparseMat(A.transpose());
    auto A_final = SparseMat(A_transpose * A);
    A_final.makeCompressed();
    auto b_final = Eigen::VectorXd(A_transpose * -b);
    solverCG.analyzePattern(A_final);
    std::cout << solverCG.info() << std::endl;
    solverCG.compute(A_final);
    std::cout << solverCG.info() << std::endl;
    //solverCG.setTolerance(.000001);
    solverCG.setMaxIterations(100000);
    Eigen::VectorXd xCG = solverCG.solve(b_final);
    //Eigen::VectorXd x = solver.solveWithGuess(b_final, Eigen::VectorXcd::Ones(num_faces*2));
    std::cout << solverCG.info() << std::endl;
    std::cout << "iters: " << solverCG.iterations() << std::endl;
    std::cout << "error: " << solverCG.error() << std::endl;
    */

    /*
    // for checking rank
    Eigen::SparseQR<SparseMat, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;
    std::cout << "size: " << A.cols() << std::endl;
    std::cout << "rank: " << solver.rank() << std::endl;
    */

    Eigen::SparseLU<SparseMat> solver(A);
    //std::cout << solver.lastErrorMessage() << std::endl;
    solver.analyzePattern(A);
    //std::cout << solver.lastErrorMessage() << std::endl;
    solver.factorize(A);
    //std::cout << "Determinant of matrix: " << solver.determinant() << std::endl;
    //std::cout << "Log of determinant of of matrix: " << solver.logAbsDeterminant() << std::endl;
    //std::cout << solver.lastErrorMessage() << std::endl;
    Eigen::VectorXd x = solver.solve(-b); // this is -b rather than b because the solver solves Ax = b, but b was set up to solve Ax + b = 0
    //std::cout << solver.info() << std::endl;
    //std::cout << solver.lastErrorMessage() << std::endl;

    // use the solution to initialize u and v for each face
    initializeDirectionFieldFromABVector(mesh, x);
}

void DirectionFieldInitializer::printSparseMatrix(const SparseMat mat) {
    for (int row = 0; row < mat.rows(); row++) {
        for (int col = 0; col < mat.cols(); col++) {
            std::cout << mat.coeff(row, col) << " ";
        }
        std::cout << std::endl;
    }
}
