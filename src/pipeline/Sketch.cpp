#include "Sketch.h"
#define NANOSVG_IMPLEMENTATION
#define NANOSVG_CPLUSPLUS
#include "lib/nanosvg/src/nanosvg.h"
#include "Mesh.h"
#include <iostream>

Eigen::Vector3f Sketch::StrokePoint::coords3d() {
    assert(triangle != nullptr);
    return triangle->vertices[0]->coords3d() * barycentric_coordinates[0] +
           triangle->vertices[1]->coords3d() * barycentric_coordinates[1] +
           triangle->vertices[2]->coords3d() * barycentric_coordinates[2];
}

float crossp(Eigen::Vector2f p1, Eigen::Vector2f p2, Eigen::Vector2f p3)
{
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
}

bool pointInTriangle(Eigen::Vector2f pt, Eigen::Vector2f v1, Eigen::Vector2f v2, Eigen::Vector2f v3)
{
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = crossp(pt, v1, v2);
    d2 = crossp(pt, v2, v3);
    d3 = crossp(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

bool lineInTriangle (Eigen::Vector2f p1, Eigen::Vector2f p2, Eigen::Vector2f v1, Eigen::Vector2f v2, Eigen::Vector2f v3) {
    Eigen::Vector3f x1(p1.x(),p1.y(),1);
    Eigen::Vector3f x2(p2.x(),p2.y(),1);
    Eigen::Matrix3f A;
    A << v1.x(), v2.x(), v3.x(),
         v1.y(), v2.y(), v3.y(),
         1, 1, 1;

    Eigen::Matrix3f inv = A.inverse();
    Eigen::Vector3f p1_barycentric = inv * x1;
    Eigen::Vector3f p2_barycentric = inv * x2;


    const float ZERO_MINUS_EP = 0 - .001;
    const float ONE_PLUS_EP = 1 + .001;
    float coeff1 = p1_barycentric.x() / (p1_barycentric.x() - p2_barycentric.x());
    Eigen::Vector3f test1 = (1 - coeff1) * p1_barycentric + coeff1 * p2_barycentric;
    if (ZERO_MINUS_EP <= coeff1 && coeff1 <= ONE_PLUS_EP && ZERO_MINUS_EP <= test1.x() && test1.x() <= ONE_PLUS_EP && ZERO_MINUS_EP <= test1.y()
            && test1.y() <= ONE_PLUS_EP && ZERO_MINUS_EP <= test1.z() && test1.z() <= ONE_PLUS_EP) {
        return true;
    }
    float coeff2 = p1_barycentric.y() / (p1_barycentric.y() - p2_barycentric.y());
    Eigen::Vector3f test2 = (1 - coeff2) * p1_barycentric + coeff2 * p2_barycentric;
    if (ZERO_MINUS_EP <= coeff2 && coeff2 <= ONE_PLUS_EP && ZERO_MINUS_EP <= test2.x() && test2.x() <= ONE_PLUS_EP && ZERO_MINUS_EP <= test2.y()
            && test2.y() <= ONE_PLUS_EP && ZERO_MINUS_EP <= test2.z() && test2.z() <= ONE_PLUS_EP) {
        return true;
    }
    float coeff3 = p1_barycentric.z() / (p1_barycentric.z() - p2_barycentric.z());
    Eigen::Vector3f test3 = (1 - coeff3) * p1_barycentric + coeff3 * p2_barycentric;
    if (ZERO_MINUS_EP <= coeff3 && coeff3 <= ONE_PLUS_EP && ZERO_MINUS_EP <= test3.x() && test3.x() <= ONE_PLUS_EP && ZERO_MINUS_EP <= test3.y()
            && test3.y() <= ONE_PLUS_EP && ZERO_MINUS_EP <= test3.z() && test3.z() <= ONE_PLUS_EP) {
        return true;
    }

    return false;
}

Eigen::Vector3f calcBarycentricCoordinates(std::shared_ptr<Face> f, std::shared_ptr<Sketch::StrokePoint> p) {
    Eigen::Vector3f x(p->coordinates.x(),p->coordinates.y(),1);
    Eigen::Matrix3f A;
    Eigen::Vector2f v1(f->vertices[0]->coords);
    Eigen::Vector2f v2(f->vertices[1]->coords);
    Eigen::Vector2f v3(f->vertices[2]->coords);
    A << v1.x(), v2.x(), v3.x(),
         v1.y(), v2.y(), v3.y(),
         1, 1, 1;
    return A.inverse() * x;
}

Sketch::Sketch(const SketchData &data) {
    width = data.width / SKETCH_SCALE;
    height = data.height / SKETCH_SCALE;
    diagonal_of_bounding_box = std::sqrt(width * width + height * height);

    for (int line_idx = 0; line_idx < data.boundary.size(); line_idx++) {
        std::vector<Eigen::Vector2f> segment;
        const std::vector<Eigen::Vector2f> &line = data.boundary[line_idx];
        for (int i = 0; i < line.size(); i++) {
            // add a boundary point and make sure it is not the same as the previous point
            if (i == 0 || (i > 0 && (line[i - 1] - line[i]).norm() > 0)) {
                segment.push_back(line[i]);
            }
        }
        boundary_strokes.push_back(std::move(segment));
    }

    int strokeId = 0;
    for (int line_idx = 0; line_idx < data.convex.size(); line_idx++) {
        Stroke segment = addBendingStrokeFromScribble(data.convex[line_idx], strokeId);
        convex_strokes.push_back(std::move(segment));
        strokeId++;
    }

    for (int line_idx = 0; line_idx < data.concave.size(); line_idx++) {
        Stroke segment = addBendingStrokeFromScribble(data.concave[line_idx], strokeId);
        concave_strokes.push_back(std::move(segment));
        strokeId++;
    }
}

Sketch::Stroke Sketch::addBendingStrokeFromScribble(const std::vector<Eigen::Vector2f> &s, int id) {
    Stroke stroke;
    stroke.id = id;

    for (int i = INDEX_OF_STROKE_TO_START; i < s.size() - INDEX_OF_STROKE_TO_START; i++) {

        // make sure this is not a duplicate of the previous point
        if (i > 0) {
            if ((s[i] - s[i-1]).norm() == 0) {
                continue;
            }
        }

        // make sure that we do not have duplicate of previous previous point
        if (i > 1) {
            if ((s[i] - s[i-2]).norm() == 0) {
                continue;
            }
        }

        std::shared_ptr<StrokePoint> point = std::make_shared<StrokePoint>();
        point->coordinates = s[i];
        point->strokeId = id;

        Eigen::Vector2f backward_point(0,0);
        int backward_points_averaged = 0;
        Eigen::Vector2f forward_point(0,0);
        int forward_points_averaged = 0;
        int offset = 1;
        while (offset <= RADIUS_OF_POINTS_TO_AVERAGE) {
            if (i - offset >= 0) { backward_point += s[i-offset]; backward_points_averaged++; }
            if (i + offset < s.size()) { forward_point += s[i+offset]; forward_points_averaged++; }
            offset++;
        }
        backward_point = backward_point / backward_points_averaged;
        forward_point = forward_point / forward_points_averaged;

        Eigen::Vector2f backward_dir(0,0);
        float backward_dist = 1;
        Eigen::Vector2f forward_dir(0,0);
        float forward_dist = 1;
        if (i == 0) {
            Eigen::Vector2f forward = forward_point - s[i];
            forward_dir = forward.normalized();
            forward_dist = forward.norm();
        } else if (i == s.size() - 1) {
            Eigen::Vector2f backward = s[i] - backward_point;
            backward_dir = backward.normalized();
            backward_dist = backward.norm();
        } else {
            Eigen::Vector2f forward = forward_point - s[i];
            forward_dir = forward.normalized();
            forward_dist = forward.norm();
            Eigen::Vector2f backward = s[i] - backward_point;
            backward_dir = backward.normalized();
            backward_dist = backward.norm();
        }

        if (backward_dist == 0 && forward_dist > 0) { point->tangent_dir = forward_dir.normalized(); }
        else if (forward_dist == 0 && backward_dist > 0) { point->tangent_dir = backward_dir.normalized(); }
        else { point->tangent_dir = (backward_dir / backward_dist + forward_dir / forward_dist).normalized(); }

        assert(point->tangent_dir.norm() > 0);
        assert(!std::isnan(point->tangent_dir.x()) && !std::isnan(point->tangent_dir.y()));

        stroke.points.push_back(point);
    }

    return stroke;
}

void Sketch::mapIntersectedFacesToStrokes(Mesh &mesh) {
    mapIntersectedFacesToStrokesHelper(mesh, concave_strokes);
    mapIntersectedFacesToStrokesHelper(mesh, convex_strokes);
}

void Sketch::mapIntersectedFacesToStrokesHelper(Mesh &mesh, std::vector<Stroke> &strokes) {
    for (int stroke_idx = 0; stroke_idx < strokes.size(); stroke_idx++) {
        int num_points = strokes[stroke_idx].points.size();
        for (int i = 0; i < num_points; i++) {
            std::shared_ptr<StrokePoint> point = strokes[stroke_idx].points[i];
            auto func = [&] (std::shared_ptr<Face> f) {
                // do a distance check so that time is not wasted on far-away triangles
                if ((f->circumcenter - point->coordinates).norm() < getTriangleAndStrokePointDistanceCheck()) {
                    // check if the StrokePoint is inside the triangle, or if one the previous
                    // line segment crosses the triangle.

                    // check if point is in triangle
                    Eigen::Vector2f v1(f->vertices[0]->coords);
                    Eigen::Vector2f v2(f->vertices[1]->coords);
                    Eigen::Vector2f v3(f->vertices[2]->coords);

                    if (pointInTriangle(point->coordinates, v1, v2, v3)) {
                        point->triangle = f;
                        point->barycentric_coordinates = calcBarycentricCoordinates(f, point);
                        if (face2strokepoints.find(f) != face2strokepoints.end()) {
                            face2strokepoints[f].push_back(point);
                        } else { face2strokepoints[f] = std::vector{point}; };
                    } else {
                        if (i > 0) {
                        // check if previous line segment crossed triangle. There is never any need to check the next line segment ahead of this point.
                            std::shared_ptr<StrokePoint> prev_point = strokes[stroke_idx].points[i-1];
                            if (lineInTriangle(point->coordinates, prev_point->coordinates, v1, v2, v3)) {
                                LineIntersection intersection;
                                intersection.points = std::make_pair(point, prev_point);
                                intersection.dir = (point->coordinates - prev_point->coordinates).normalized();
                                assert(!std::isnan(intersection.dir.x()) && !std::isnan(intersection.dir.y()));
                                assert(!std::isinf(intersection.dir.x()) && !std::isinf(intersection.dir.y()));
                                if (face2lines.find(f) != face2lines.end()) {
                                    face2lines[f].push_back(intersection);
                                } else {
                                    face2lines[f] = std::vector{intersection};
                                }
                            }
                        }
                    }
                }
            };

            mesh.forEachTriangle(func);
            if (point->triangle == nullptr) {
                std::cout << "WARNING: there is a bending stroke point that is not inside of a triangle!" << std::endl;
                // TODO: do something better about this
            }
        }
    }
}
