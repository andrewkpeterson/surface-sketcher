#include "HeightFieldSolver.h"
#include "autodiff/forward.hpp"
#include "alglib/optimization.h"
#include <stdlib.h>
#include <stdio.h>
#include "ceres/ceres.h"
#include "src/utils/OBJWriter.h"

void HeightFieldSolver::solveForHeightField(Mesh &mesh, Sketch &sketch) {

    /*
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        f->vertices[0]->height = std::sin(f->vertices[0]->coords.y() * 2);
        f->vertices[1]->height = std::sin(f->vertices[1]->coords.y() * 2);
        f->vertices[2]->height = std::sin(f->vertices[2]->coords.y() * 2);
    });
    estimateCurvatureValues(mesh, sketch);
    */

    for (int i = 0; i < NUM_ITERATIONS; i++) {
        if (i == 0) {
            initializeCurvatureValues(sketch);
        } else {
            estimateCurvatureValues(mesh, sketch);
        }
        minimizeELambda(mesh, sketch);
        optimizeHeightField(mesh, sketch);
    }

}

void HeightFieldSolver::initializeCurvatureValues(Sketch &sketch) {

    auto &concave_strokes = sketch.getConcaveStrokes();
    auto &convex_strokes = sketch.getConvexStrokes();

    // Convex sections (i.e. "Hills") have negative normal curvature because the unit tangent vector
    // pulls in the OPPOSITE direction of the surface normal, which points upward. Concave sections
    // (i.e. "Valleys") have positive normal curvature because the unit tangent vector pulls into the
    // SAME direction as the surface normal.
    for (int stroke_idx = 0; stroke_idx < concave_strokes.size(); stroke_idx++) {
        auto &stroke = concave_strokes[stroke_idx];
        int num_points = sketch.getLengthOfStrokeInPoints(stroke);
        for (int i = 0; i < num_points; i++) {
            sketch.getStrokePointByFlattenedIndex(stroke, i)->curvature_value = INITIAL_CURVATURE_MAGNITUDE;
        }
    }

    for (int stroke_idx = 0; stroke_idx < convex_strokes.size(); stroke_idx++) {
        auto &stroke = convex_strokes[stroke_idx];
        int num_points = sketch.getLengthOfStrokeInPoints(stroke);
        for (int i = 0; i < num_points; i++) {
            sketch.getStrokePointByFlattenedIndex(stroke, i)->curvature_value = -INITIAL_CURVATURE_MAGNITUDE;
        }
    }

}

void HeightFieldSolver::estimateCurvatureValues(Mesh &mesh, Sketch &sketch) {
    estimateCurvatureValuesHelper(mesh, sketch, sketch.getConvexStrokes(), true);
    estimateCurvatureValuesHelper(mesh, sketch, sketch.getConcaveStrokes(), false);

}

struct F1 {
  F1(double x, double y) : p1(x), p2(y) {}
  template <typename T>
  bool operator()(const T* const c1, const T* const c2, const T* const r,  T* residual) const {
    residual[0] = T(sqrt((T(p1) - c1[0]) * (T(p1) - c1[0]) + (T(p2) - c2[0]) * (T(p2) - c2[0]))) - r[0];
    //std::cout << residual[0] << std::endl;
    return true;
  }

private:
  double p1;
  double p2;
};

struct F2 {
  template <typename T>
  bool operator()(const T* const r, T* residual) const {
    residual[0] = T(std::sqrt(HeightFieldSolver::PRINCIPLE_CURVATURE_MU)) * r[0];
    //residual[0] = T(HeightFieldSolver::PRINCIPLE_CURVATURE_MU) * r[0];
    return true;
  }
};

void HeightFieldSolver::estimateCurvatureValuesHelper(Mesh &mesh, Sketch &sketch,
                                                      std::vector<Sketch::Stroke> &strokes, bool convex) {
    for (int stroke_idx = 0; stroke_idx < strokes.size(); stroke_idx++) {
        for (int seg_idx = 0; seg_idx < strokes[stroke_idx].segments.size(); seg_idx++) {
            auto &segment = strokes[stroke_idx].segments[seg_idx].seg;
            std::vector<Eigen::Vector2f> points;
            int center_idx = segment.size() / 2;
            assert(segment[center_idx]->triangle != nullptr);

            // Here, we are creating the normal section plane that defines the curve on the surface whose
            // normal curvature we are trying to calculate. Here, we lift the stroke points to the surface
            // and then fit a circle to the planar curve on the normal section plane in order to estimate the
            // curvature of the surface.
            Eigen::Vector3f origin = segment[center_idx]->coords3d();
            Eigen::Vector3f n = segment[center_idx]->triangle->normal().normalized();
            Eigen::Vector3f t = Eigen::Vector3f(segment[center_idx]->tangent_dir[0], segment[center_idx]->tangent_dir[1], 0);
            Eigen::Vector3f plane_normal = t.cross(n).normalized();
            Eigen::Vector3f other_basis = n.cross(plane_normal).normalized();
            for (int p_idx = 0; p_idx < segment.size(); p_idx++) {
                // NOTE: This will work best when the strokes are very dense with points
                // NOTE: This might work better if the curvature was estimated at every
                // stroke point, rather than on each stroke point segment
                Eigen::Vector3f projected_point = segment[p_idx]->coords3d() - (segment[p_idx]->coords3d() - origin).dot(plane_normal) * plane_normal;
                Eigen::Vector2f point_on_plane = Eigen::Vector2f((projected_point - origin).dot(n), (projected_point - origin).dot(other_basis));
                points.push_back(point_on_plane);
            }

            using ceres::AutoDiffCostFunction;
            using ceres::CostFunction;
            using ceres::Problem;
            using ceres::Solve;
            using ceres::Solver;

            // We need to initialize the center of the circle according to the known convexity of the bending line.
            // "Hills" in the surface are due to convex bending lines, and the surface normal always points upward.
            // This means that the center of the circle we are creating will be BELOW the surface (so we initialize
            // the center of the circle at (-.1,-.1) on the projection plane spanned by the lifted curve's tangent
            // vector and the surface normal). "Valleys" in the surface are due to conave bending lines, which means
            // than the center of the circle should be ABOVE the projection plane. This is why the center of the
            // circle is initialized at (.1,.1)
            // Remember that the origin of the projection plane is the center stroke point of the line segment, so we
            // don't actually need to take the height of the surface into account when initializing the center.
            double r = 10.0;
            double c1 = convex ? -.1 : .1;
            double c2 = convex ? -.1 : .1;

            Problem problem;

            for (int i = 0; i < segment.size(); i++) {
                CostFunction* func1 = new AutoDiffCostFunction<F1, 1, 1, 1, 1>(new F1(points[i].x(), points[i].y()));
                CostFunction* func2 = new AutoDiffCostFunction<F2, 1, 1>(new F2);
                problem.AddResidualBlock(func1, nullptr, &c1, &c2, &r);
                problem.AddResidualBlock(func2, nullptr, &r);
            }

            Solver::Options options;
            options.max_num_iterations = 200;
            options.function_tolerance = .000001;
            //options.minimizer_progress_to_stdout = true;
            Solver::Summary summary;
            Solve(options, &problem, &summary);
            std::cout << summary.BriefReport() << "\n";

            // Update curvature values for segment and points. Remember that we are measuring normal curvature!
            // Convex sections (i.e. "Hills") have negative normal curvature because the unit tangent vector
            // pulls in the opposite direction of the surface normal, which points upward. Concave sections
            // (i.e. "Valleys") have positive normal curvature because the unit tangent vector pull into the
            // same direction as the surface normal.
            for (int i = 0; i < segment.size(); i++) {
                segment[i]->curvature_value = convex ? -(1.0 / r) : (1.0 / r);
            }
        }
    }
}

void HeightFieldSolver::minimizeELambda(Mesh &mesh, const Sketch &sketch) {

    int num_faces = mesh.getNumTriangles();

    // get A and b to calculate new u vectors
    std::vector<Triplet> A_coefficients;
    SparseMat A(num_faces*2, num_faces*2);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(num_faces*2);
    std::map<std::pair<int,int>, double> A_map;

    addCoefficientsForELambda(mesh, sketch, A_map, true, true);
    addCoefficientsForELambda(mesh, sketch, A_map, true, false);
    addCoefficientsForELambda(mesh, sketch, A_map, false, true);
    addCoefficientsForELambda(mesh, sketch, A_map, false, false);
    addConstraintsForELambda(mesh, sketch, A_map, b);

    for (auto it = A_map.begin(); it != A_map.end(); it++) {
        A_coefficients.push_back(Triplet(it->first.first, it->first.second, it->second));
    }
    A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());

    // solve for u vectors
    A.makeCompressed();
    Eigen::SparseLU<SparseMat> solver(A);
    solver.analyzePattern(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;
    solver.factorize(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;
    Eigen::VectorXd x = solver.solve(-b); // this is -b rather than bu because the solver solves Ax = bu, but b was set up to solve Ax + bu = 0
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;

    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        f->lambda_u = x(2*f->index);
        f->lambda_v = x(2*f->index + 1);
    });
}

void HeightFieldSolver::addCoefficientsForELambda(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>,double> &m, bool constraining_lambda_u, bool in_dir_u) {

    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        Eigen::Vector2f p = f->circumcenter;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 3);
        for (int i = 0; i < f->neighbors.size(); i++) {
            A(0,i) = f->neighbors[i]->circumcenter.x() - p.x();
            A(1,i) = f->neighbors[i]->circumcenter.y() - p.y();
        }

        Eigen::Vector2f q = in_dir_u ? f->u : f->v;
        Eigen::MatrixXd A_pinv = A.completeOrthogonalDecomposition().pseudoInverse();

        float k = q.x() * (A_pinv(0,0) + A_pinv(1,0) + A_pinv(2,0)) + q.y() * (A_pinv(0,1) + A_pinv(1,1) + A_pinv(2,1));
        float k1 = -(q.x() * A_pinv(0,0) + q.y() * A_pinv(0,1));
        float k2 = -(q.x() * A_pinv(1,0) + q.y() * A_pinv(1,1));
        float k3 = -(q.x() * A_pinv(2,0) + q.y() * A_pinv(2,1));

        assert((k1 == 0 && f->neighbors.size() < 1) || f->neighbors.size() >= 1);
        assert((k2 == 0 && f->neighbors.size() < 2) || f->neighbors.size() >= 2);
        assert((k3 == 0 && f->neighbors.size() < 3) || f->neighbors.size() >= 3);

        bool need_smoothness_beta = (constraining_lambda_u && in_dir_u) || (!constraining_lambda_u && !in_dir_u);

        float mult = 2 * f->area / mesh.getTotalArea() * (need_smoothness_beta ? CURVATURE_MAGNITUDE_SMOOTHNESS_BETA : 1);

        {
            int idx = 2*f->index + (constraining_lambda_u ? 0 : 1);
            int idx1 = 2*(f->neighbors.size() > 0 ? f->neighbors[0]->index : 0) + (constraining_lambda_u ? 0 : 1);
            int idx2 = 2*(f->neighbors.size() > 1 ? f->neighbors[1]->index : 0) + (constraining_lambda_u ? 0 : 1);
            int idx3 = 2*(f->neighbors.size() > 2 ? f->neighbors[2]->index : 0) + (constraining_lambda_u ? 0 : 1);

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

void HeightFieldSolver::addConstraintsForELambda(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>,double> &m, Eigen::VectorXd &b) {

    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        std::set<int> strokes;
        if (sketch.checkStrokePoints(f)) {
            for (int i = 0; i < sketch.getStrokePoints(f).size(); i++) {
                if (strokes.find(sketch.getStrokePoints(f)[i]->strokeId) == strokes.end()) {
                    Eigen::Vector2f dir = sketch.getStrokePoints(f)[i]->tangent_dir.normalized();
                    addConstraintsForELambdaHelper(mesh, sketch, m, f, b, dir, sketch.getStrokePoints(f)[i]->curvature_value);
                }
                strokes.insert(sketch.getStrokePoints(f)[i]->strokeId);
            }
        }

        if (sketch.checkStrokeLines(f)) {
            for (int i = 0; i < sketch.getStrokeLines(f).size(); i++) {
                if (strokes.find(sketch.getStrokeLines(f)[i].points.first->strokeId) == strokes.end()) {
                    Eigen::Vector2f dir = sketch.getStrokeLines(f)[i].dir.normalized();
                    float curvature_value = (sketch.getStrokeLines(f)[i].points.first->curvature_value +
                                             sketch.getStrokeLines(f)[i].points.second->curvature_value) / 2.0;
                    addConstraintsForELambdaHelper(mesh, sketch, m, f, b, dir, curvature_value);
                }
                strokes.insert(sketch.getStrokeLines(f)[i].points.first->strokeId);
            }
        }
    });

}

void HeightFieldSolver::addConstraintsForELambdaHelper(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int, int>, double> &m, std::shared_ptr<Face> f,
                                                       Eigen::VectorXd &b, Eigen::Vector2f d, float curvature_value) {

    Eigen::Vector2f best_u = (f->u.dot(d)) > (-f->u.dot(d)) ? f->u: -f->u;
    Eigen::Vector2f best_v = (f->v.dot(d)) > (-f->v.dot(d)) ? f->v : -f->v;
    bool forV = best_v.dot(d) > best_u.dot(d);

    int idx = forV ? 2*f->index + 1 : 2*f->index;

    float mult = 2 * f->area / mesh.getTotalArea() * CURVATURE_MAGNITUDE_CONSTRAINT_WEIGHT;
    std::pair<int, int> x_idx(idx, idx);
    addToSparseMap(x_idx, mult, m);
    b(idx) += -mult * curvature_value;

}

void HeightFieldSolver::optimizeHeightField(Mesh &mesh, const Sketch &sketch) {

    int num_vertices = mesh.getNumVertices();

    // get A and b to calculate new u vectors
    std::vector<Triplet> A_coefficients;
    SparseMat A(num_vertices, num_vertices);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(num_vertices);
    std::map<std::pair<int,int>, double> A_map;

    optimizeHeightFieldHelper(mesh, sketch, A_map, b);

    for (auto it = A_map.begin(); it != A_map.end(); it++) {
        A_coefficients.push_back(Triplet(it->first.first, it->first.second, it->second));
    }
    A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());

    // solve for u vectors
    A.makeCompressed();
    Eigen::SparseLU<SparseMat> solver(A);
    solver.analyzePattern(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;
    solver.factorize(A);
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;
    Eigen::VectorXd x = solver.solve(-b); // this is -b rather than bu because the solver solves Ax = bu, but b was set up to solve Ax + bu = 0
    std::cout << solver.lastErrorMessage() << std::endl;
    std::cout << solver.info() << std::endl;

    mesh.forEachVertex([&](std::shared_ptr<Vertex> v) {
        v->height = x(v->index);
    });

}

void HeightFieldSolver::optimizeHeightFieldHelper(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b) {

    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {

        if (f->neighbors.size() < 3) { return; }

        // add coefficients for E_match
        {
            Eigen::Vector2f p = f->circumcenter;
            Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
            for (int i = 0; i < f->neighbors.size(); i++) {
                Eigen::Vector2f n = f->neighbors[i]->circumcenter;
                A(i,0) = std::pow(n.x() - p.x(), 2);
                A(i,1) = (n.x() - p.x()) * (n.y() - p.y());
                A(i,2) = std::pow(n.y() - p.y(), 2);
            }

            auto A_inv = A.inverse();
            auto G_inv = calculateGInverse(f);

            HeightFieldSolver::Values vs;

            vs.n3 = f->normal()(2);

            vs.G11 = G_inv(0,0); vs.G12 = G_inv(0,1);
            vs.G21 = G_inv(1,0); vs.G22 = G_inv(1,1);

            vs.a11 = A_inv(0,0); vs.a12 = A_inv(0,1); vs.a13 = A_inv(0,2);
            vs.a21 = A_inv(1,0); vs.a22 = A_inv(1,1); vs.a23 = A_inv(1,2);
            vs.a31 = A_inv(2,0); vs.a32 = A_inv(2,1); vs.a33 = A_inv(2,2);

            vs.faces = {f.get(), f->neighbors[0], f->neighbors[1], f->neighbors[2]};

            computeAndAddCoefficientsForEMatch(mesh, sketch, m, b, f, true, true, vs);
            computeAndAddCoefficientsForEMatch(mesh, sketch, m, b, f, true, false, vs);
            computeAndAddCoefficientsForEMatch(mesh, sketch, m, b, f, false, true, vs);
            computeAndAddCoefficientsForEMatch(mesh, sketch, m, b, f, false, false, vs);
        }

        // add coefficients to E_bdry (i.e. boundary conditions)
        {
            addCoefficientsForEBoundary(m, b);
        }

    });

}

void HeightFieldSolver::computeAndAddCoefficientsForEMatch(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b,
                                                           std::shared_ptr<Face> f, bool u_direction, bool x_coordinate, HeightFieldSolver::Values &vs) {

    float K_1 = calcK(f.get(), f->vertices[0], 0, vs, u_direction, x_coordinate);
    float K_2 = calcK(f.get(), f->vertices[1], 0, vs, u_direction, x_coordinate);
    float K_3 = calcK(f.get(), f->vertices[2], 0, vs, u_direction, x_coordinate);

    float K1_1 = calcK(f->neighbors[0], f->neighbors[0]->vertices[0], 1, vs, u_direction, x_coordinate);
    float K1_2 = calcK(f->neighbors[0], f->neighbors[0]->vertices[1], 1, vs, u_direction, x_coordinate);
    float K1_3 = calcK(f->neighbors[0], f->neighbors[0]->vertices[2], 1, vs, u_direction, x_coordinate);

    float K2_1 = calcK(f->neighbors[1], f->neighbors[1]->vertices[0], 2, vs, u_direction, x_coordinate);
    float K2_2 = calcK(f->neighbors[1], f->neighbors[1]->vertices[1], 2, vs, u_direction, x_coordinate);
    float K2_3 = calcK(f->neighbors[1], f->neighbors[1]->vertices[2], 2, vs, u_direction, x_coordinate);

    float K3_1 = calcK(f->neighbors[2], f->neighbors[2]->vertices[0], 3, vs, u_direction, x_coordinate);
    float K3_2 = calcK(f->neighbors[2], f->neighbors[2]->vertices[1], 3, vs, u_direction, x_coordinate);
    float K3_3 = calcK(f->neighbors[2], f->neighbors[2]->vertices[2], 3, vs, u_direction, x_coordinate);

    std::vector<float> coefficients = {K_1, K_2, K_3, K1_1, K1_2, K1_3, K2_1, K2_2, K2_3, K3_1, K3_2, K3_3};
    std::vector<std::shared_ptr<Vertex>> vertices = {f->vertices[0], f->vertices[1], f->vertices[2],
                                                     f->neighbors[0]->vertices[0], f->neighbors[0]->vertices[1], f->neighbors[0]->vertices[2],
                                                     f->neighbors[1]->vertices[0], f->neighbors[1]->vertices[1], f->neighbors[1]->vertices[2],
                                                     f->neighbors[2]->vertices[0], f->neighbors[2]->vertices[1], f->neighbors[2]->vertices[2]};

    float constraint_lambda = u_direction ? f->lambda_u : f->lambda_v;
    Eigen::Vector2f constraint_direction = u_direction ? f->u : f->v;
    float constraint_coordinate = x_coordinate ? constraint_direction.x() : constraint_direction.y();
    float constraint = -constraint_lambda * constraint_coordinate;

    addCoefficientsToMapAndVectorForEMatch(coefficients, vertices, m, b, constraint);

}

Eigen::Matrix2f HeightFieldSolver::calculateGInverse(std::shared_ptr<Face> f) {
    std::shared_ptr<Vertex> i = f->vertices[0];
    std::shared_ptr<Vertex> j = f->vertices[1];
    std::shared_ptr<Vertex> k = f->vertices[2];

    Eigen::Matrix2f m;
    m << 0, -1,
         1, 0;
    Eigen::Matrix2f m_reverse;
    m << 0, 1,
         -1, 0;

    Eigen::Vector2f ki = i->coords - k->coords;
    Eigen::Vector2f ki_perp = (m * ki).dot(j->coords - i->coords) > 0 ? m * ki : m_reverse * ki;

    Eigen::Vector2f ij = j->coords - i->coords;
    Eigen::Vector2f ij_perp = (m * ij).dot(k->coords - j->coords) > 0 ? m * ij : m_reverse * ij;

    Eigen::Vector2f derivs = ((j->height - i->height) * ki_perp / (2.0 * f->area)) + ((k->height - i->height) * ij_perp / 2.0 * f->area);

    Eigen::Matrix2f G;
    G << (1 + derivs.x() * derivs.x()), (derivs.x() * derivs.y()),
         (derivs.x() * derivs.y()), (1 + derivs.y() * derivs.y());

    return G.inverse();
}

float HeightFieldSolver::calcK(Face* f, std::shared_ptr<Vertex> v, int face_idx, HeightFieldSolver::Values &vs,
                               bool u_direction, bool x_coordinate) {
    Eigen::Vector2f w = u_direction ? f->u : f->v;
    float coefficient = x_coordinate ? (w.x() * calcdN(f,v,face_idx,vs,1) + w.y() * calcdN(f,v,face_idx,vs,2)) :
                                       (w.x() * calcdN(f,v,face_idx,vs,3) + w.y() * calcdN(f,v,face_idx,vs,4));
    return coefficient;
}

float HeightFieldSolver::calcdN(Face* f, std::shared_ptr<Vertex> v, int face_idx, HeightFieldSolver::Values &vs, int N_idx) {
    switch (N_idx) {
        case 1: return (vs.G11 * calcdAlpha(f,v,face_idx,vs) + vs.G21 * calcdBeta(f,v,face_idx,vs)) * -vs.n3;
        case 2: return (vs.G12 * calcdAlpha(f,v,face_idx,vs) + vs.G22 * calcdBeta(f,v,face_idx,vs)) * -vs.n3;
        case 3: return (vs.G11 * calcdBeta(f,v,face_idx,vs) + vs.G21 * calcdGamma(f,v,face_idx,vs)) * -vs.n3;
        case 4: return (vs.G12 * calcdBeta(f,v,face_idx,vs) + vs.G22 * calcdGamma(f,v,face_idx,vs)) * -vs.n3;
    }
}

float HeightFieldSolver::calcdAlpha(Face* f, std::shared_ptr<Vertex> v, int face_idx, HeightFieldSolver::Values &vs) {
    return vs.a11 * calcdB(f,v,face_idx,1) + vs.a12 * calcdB(f,v,face_idx,2) + vs.a13 * calcdB(f,v,face_idx,3);
}

float HeightFieldSolver::calcdBeta(Face* f, std::shared_ptr<Vertex> v, int face_idx, HeightFieldSolver::Values &vs) {
    return vs.a21 * calcdB(f,v,face_idx,1) + vs.a22 * calcdB(f,v,face_idx,2) + vs.a23 * calcdB(f,v,face_idx,3);
}

float HeightFieldSolver::calcdGamma(Face* f, std::shared_ptr<Vertex> v, int face_idx, Values &vs) {
    return vs.a31 * calcdB(f,v,face_idx,1) + vs.a32 * calcdB(f,v,face_idx,2) + vs.a33 * calcdB(f,v,face_idx,3);
}

float HeightFieldSolver::calcdB(Face* f, std::shared_ptr<Vertex> v, int face_idx, int B_idx) {
    if (face_idx != 0 && face_idx != B_idx) { return 0; }
    Face *n = f->neighbors[B_idx - 1]; // B_idx == 1 + neighbor_index
    return (n->circumcenter.x() - f->circumcenter.x()) * calcdG(f,v,face_idx,true) +
           (n->circumcenter.y() - f->circumcenter.y()) * calcdG(f,v,face_idx,false);
}

float HeightFieldSolver::calcdG(Face* f, std::shared_ptr<Vertex> v, int face_idx, bool x_coordinate) {
    std::shared_ptr<Vertex> i;
    std::shared_ptr<Vertex> j;
    std::shared_ptr<Vertex> k;
    if (v == f->vertices[0]) { i = v; j = f->vertices[1]; k = f->vertices[2]; }
    if (v == f->vertices[1]) { i = v; j = f->vertices[0]; k = f->vertices[2]; }
    if (v == f->vertices[2]) { i = v; j = f->vertices[0]; k = f->vertices[1]; }

    Eigen::Matrix2f m;
    m << 0, -1,
         1, 0;
    Eigen::Matrix2f m_reverse;
    m << 0, 1,
         -1, 0;

    Eigen::Vector2f ki = i->coords - k->coords;
    Eigen::Vector2f ki_perp = (m * ki).dot(j->coords - i->coords) > 0 ? m * ki : m_reverse * ki;

    Eigen::Vector2f ij = j->coords - i->coords;
    Eigen::Vector2f ij_perp = (m * ij).dot(k->coords - j->coords) > 0 ? m * ij : m_reverse * ij;

    float C_ik = x_coordinate ? (ki_perp.x() / (2.0 * f->area)) : (ki_perp.y() / (2.0 * f->area));
    float C_ji = x_coordinate ? (ij_perp.x() / (2.0 * f->area)) : (ij_perp.y() / (2.0 * f->area));

    return (face_idx == 0) ? (C_ik + C_ji) : -(C_ik + C_ji);
}

void HeightFieldSolver::addCoefficientsToMapAndVectorForEMatch(std::vector<float> coefficients, std::vector<std::shared_ptr<Vertex>> vertices,
                                                       std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b, float constraint) {

    for (int vertex_idx = 0; vertex_idx < coefficients.size(); vertex_idx++) {
        auto v = vertices[vertex_idx]; // the vertex that the with respect to which the derivative was taken
        float kv = coefficients[vertex_idx];
        b(v->index) += constraint * kv;
        for (int neighbor_idx = 0; neighbor_idx < coefficients.size(); neighbor_idx++) {
            auto n = vertices[neighbor_idx];
            float kn = coefficients[neighbor_idx];
            std::pair<int,int> p(v->index, n->index);
            addToSparseMap(p, 2 * kv * kn, m);
        }
    }

}

void HeightFieldSolver::addCoefficientsForEBoundary(Mesh &mesh, const Sketch &sketch, std::map<std::pair<int,int>, double> &m, Eigen::VectorXd &b) {
    mesh.forEachBoundaryVertex([&](std::shared_ptr<const Vertex> v) {

    });
}
