#include "HeightFieldSolver.h"
#include "autodiff/forward.hpp"
#include "alglib/optimization.h"
#include <stdlib.h>
#include <stdio.h>
#include "ceres/ceres.h"

void HeightFieldSolver::solveForHeightField(Mesh &mesh, Sketch &sketch) {

    mesh.forEachTriangle([&](std::shared_ptr<Face> f) {
        f->vertices[0]->height = std::sin(f->vertices[0]->coords.y() * 2);
        f->vertices[1]->height = std::sin(f->vertices[1]->coords.y() * 2);
        f->vertices[2]->height = std::sin(f->vertices[2]->coords.y() * 2);
    });
    estimateCurvatureValues(mesh, sketch);

    /*
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        if (i == 0) {
            initializeCurvatureValues(sketch);
        } else {
            estimateCurvatureValues(mesh, sketch);
        }
        minimizeELambda(mesh, sketch);
        minimizeEMatch(mesh, sketch);
    }
    */

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
                /*
                std::cout << "coords: " << segment[p_idx]->coords3d().x() << " " << segment[p_idx]->coords3d().y() << " " << segment[p_idx]->coords3d().z() << std::endl;
                std::cout << "projected point: " << point_on_plane.x() << " " << point_on_plane.y() << std::endl;
                std::cout << "plane_normal:" << plane_normal.x() << " " << plane_normal.y() << " " << plane_normal.z() << std::endl;
                std::cout << "n:" << n.x() << " " << n.y() << " " << n.z() << std::endl;
                std::cout << "t:" << t.x() << " " << t.y() << " " << t.z() << std::endl;
                */
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



}

void HeightFieldSolver::minimizeEMatch(Mesh &mesh, const Sketch &sketch) {



}
