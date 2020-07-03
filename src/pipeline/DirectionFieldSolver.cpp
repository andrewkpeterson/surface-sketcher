#include "DirectionFieldSolver.h"
#include "Mesh.h"
#include <Eigen/Sparse>
#include <map>



void DirectionFieldSolver::initializeDirectionField(Mesh &mesh, Sketch &sketch) {

    using SparseMat = Eigen::SparseMatrix<std::complex<float>>;
    using CTriplet = Eigen::Triplet<std::complex<float>>;

    int num_faces = mesh.getNumTriangles();

    SparseMat A(num_faces*2, num_faces*2); // multiply faces by 2 to account for a and b values
    Eigen::VectorXcf b(num_faces*2);

    std::map<std::pair<int,int>, std::complex<float>> sparse_map;

    // add coefficients for E_smooth
    mesh.forEachConstPairOfNeighboringTriangles([&](const Face *f1, const Face *f2) {
        float efg = Mesh::calcEFGArea(f1, f2);
        if ()
    });

    // add coefficients for E_constraint
    mesh.forEachConstTriangle([&](std::shared_ptr<const Face> f) {

    });

    // add coefficients for E_ortho
    mesh.forEachConstTriangle([&](std::shared_ptr<const Face> f) {

    });

    // turn the map into triples and put the triples into the sparse matrix

    // solve the system of linear equations

    // use the solution to initialize u and v for each face

    /*
    std::vector<CTriplet> triplets;
    //triplets.reserve(estimation_of_entries);
    for (auto it = mesh.cdt.faces_begin(); it != mesh.cdt.faces_end(); it++) {
        triplets.push_back(CTriplet());
    }
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat.
    */

    // add coefficients for E_ortho


}

void DirectionFieldSolver::optimizeBendFieldEnergy(Mesh &mesh, Sketch &sketch) {

}
