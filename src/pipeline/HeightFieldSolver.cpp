#include "HeightFieldSolver.h"
#include "autodiff/forward.hpp"
#include "alglib/optimization.h"
#include <stdlib.h>
#include <stdio.h>
#include "ceres/ceres.h"

void HeightFieldSolver::solveForHeightField(Mesh &mesh, Sketch &sketch) {

    /*
    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        f->vertices[0]->height = .01 * f->vertices[0]->coords.x() * f->vertices[0]->coords.y();
        f->vertices[1]->height = .01 * f->vertices[1]->coords.x() * f->vertices[1]->coords.y();
        f->vertices[2]->height = .01 * f->vertices[2]->coords.x() * f->vertices[2]->coords.y();
    });
    */
    estimateCurvatureValues(mesh, sketch);
    /*
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        if (i == 0) {
            initializeCurvatureValues(sketch);
        } else {
            estimateCurvatureValues(mesh, sketch);
        }
        estimateCurvatureValues(mesh, sketch);
        minimizeELambda(mesh, sketch);
        minimizeEMatch(mesh, sketch);
    }
    */

}

void HeightFieldSolver::initializeCurvatureValues(Sketch &sketch) {

    auto &concave_strokes = sketch.getConcaveStrokes();
    auto &convex_strokes = sketch.getConvexStrokes();

    for (int stroke_idx = 0; stroke_idx < concave_strokes.size(); stroke_idx++) {
        for (int seg_idx = 0; seg_idx < concave_strokes[stroke_idx].size(); seg_idx++) {
            concave_strokes[stroke_idx][seg_idx].curvature_value = -INITIAL_CURVATURE_MAGNITUDE;
            for (int i = 0; i < concave_strokes[stroke_idx][seg_idx].seg.size(); i++) {
                concave_strokes[stroke_idx][seg_idx].seg[i]->curvature_value = -INITIAL_CURVATURE_MAGNITUDE;
            }
        }
    }

    for (int stroke_idx = 0; stroke_idx < convex_strokes.size(); stroke_idx++) {
        for (int seg_idx = 0; seg_idx < convex_strokes[stroke_idx].size(); seg_idx++) {
            convex_strokes[stroke_idx][seg_idx].curvature_value = INITIAL_CURVATURE_MAGNITUDE;
            for (int i = 0; i < convex_strokes[stroke_idx][seg_idx].seg.size(); i++) {
                convex_strokes[stroke_idx][seg_idx].seg[i]->curvature_value = INITIAL_CURVATURE_MAGNITUDE;
            }
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
    return true;
  }
};

void HeightFieldSolver::estimateCurvatureValuesHelper(Mesh &mesh, Sketch &sketch,
                                                      std::vector<std::vector<Sketch::CurvatureStrokeSegment>> &strokes, bool convex) {
    for (int stroke_idx = 0; stroke_idx < strokes.size(); stroke_idx++) {
        for (int seg_idx = 0; seg_idx < strokes[stroke_idx].size(); seg_idx++) {
            auto &segment = strokes[stroke_idx][seg_idx].seg;
            std::vector<Eigen::Vector2f> points;
            int center_idx = segment.size() / 2;
            assert(segment[center_idx]->triangle != nullptr);
            Eigen::Vector3f origin = segment[center_idx]->coords3d();
            Eigen::Vector3f n = segment[center_idx]->triangle->normal().normalized();
            Eigen::Vector3f t = Eigen::Vector3f(segment[center_idx]->tangent_dir[0], segment[center_idx]->tangent_dir[1], 0);
            Eigen::Vector3f plane_normal = t.cross(n).normalized();
            Eigen::Vector3f other_basis = n.cross(plane_normal).normalized();
            for (int p_idx = 0; p_idx < segment.size(); p_idx++) {
                // lift the stroke to the surface
                // this probably works best when the strokes are very dense with points
                Eigen::Vector3f projected_point = segment[p_idx]->coords3d() - (segment[p_idx]->coords3d() - origin).dot(plane_normal) * plane_normal;
                Eigen::Vector2f point_on_plane = Eigen::Vector2f((projected_point - origin).dot(n), (projected_point - origin).dot(other_basis));
                points.push_back(point_on_plane);
            }

            using ceres::AutoDiffCostFunction;
            using ceres::CostFunction;
            using ceres::Problem;
            using ceres::Solve;
            using ceres::Solver;

            double r = 1.0;
            double c1 = -.1;
            double c2 = -.1;

            Problem problem;

            for (int i = 0; i < segment.size(); i++) {
                CostFunction* func1 = new AutoDiffCostFunction<F1, 1, 1, 1, 1>(new F1(points[i].x(), points[i].y()));
                CostFunction* func2 = new AutoDiffCostFunction<F2, 1, 1>(new F2);
                problem.AddResidualBlock(func1, nullptr, &c1, &c2, &r);
                problem.AddResidualBlock(func2, nullptr, &r);
            }

            Solver::Options options;
            options.max_num_iterations = 200;
            options.function_tolerance = .000000001;
            //options.minimizer_progress_to_stdout = true;
            Solver::Summary summary;
            Solve(options, &problem, &summary);
            std::cout << summary.BriefReport() << "\n";

            // update curvature values for segment and points
            strokes[stroke_idx][seg_idx].curvature_value = convex ? (1.0 / r) : -(1.0 / r);
            for (int i = 0; i < segment.size(); i++) {
                segment[i]->curvature_value = convex ? (1.0 / r) : -(1.0 / r);
            }
            std::cout << "******************** " << strokes[stroke_idx][seg_idx].curvature_value << " *********************" << std::endl;
        }
    }
}

void HeightFieldSolver::minimizeELambda(Mesh &mesh, const Sketch &sketch) {



}

void HeightFieldSolver::minimizeEMatch(Mesh &mesh, const Sketch &sketch) {



}
