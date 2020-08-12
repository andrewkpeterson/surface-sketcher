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

std::shared_ptr<Sketch::StrokePoint> Sketch::getStrokePointByFlattenedIndex(const Stroke &stroke, int idx) {
    assert(idx < getLengthOfStrokeInPoints(stroke));
    for (int seg_idx = 0; seg_idx < stroke.segments.size(); seg_idx++) {
        for (int point_idx = 0; point_idx < stroke.segments[seg_idx].seg.size(); point_idx++) {
            if (idx <= 0) {
                return stroke.segments[seg_idx].seg[point_idx];
            } else {
                idx--;
            }
        }
    }
    assert(false);
}

int Sketch::getLengthOfStrokeInPoints(const Stroke &stroke) {
    int length = 0;
    for (int i = 0; i < stroke.segments.size(); i++) {
        length += stroke.segments[i].seg.size();
    }
    return length;
}

Sketch::Sketch(std::string svg_file) {
    NSVGimage *image;
    image = nsvgParseFromFile(svg_file.c_str(), "px", SVG_SIZE_PARAM);

    width = image->width;
    height = image->height;
    diagonal_of_bounding_box = std::sqrt(width * width + height * height);

    NSVGshape *shape;
    for (shape = image->shapes; shape != NULL; shape = shape->next) {
        if (shape->stroke.color == 0xFF000000) {
            NSVGpath *path = shape->paths;
            std::vector<Eigen::Vector2f> line;
            for (int i = 0; i < path->npts - 1; i += 2) {
                // add a boundary point and make sure it is not the same as the previous point
                if (i == 0 || (i > 0 && (Eigen::Vector2f(path->pts[i], path->pts[i+1]) - Eigen::Vector2f(path->pts[i-2], path->pts[i-1])).norm() > 0)) {
                    line.push_back(Eigen::Vector2f(path->pts[i], path->pts[i+1]));
                }
            }
            boundary_strokes.push_back(line);
        } else if (shape->stroke.color == 0xFF0000FF) { // red convex bending line
            NSVGpath *path = shape->paths;
            Stroke stroke = addBendingStrokeFromSVG(path);
            convex_strokes.push_back(stroke);
        } else if (shape->stroke.color == 0xFFFF0000) { // blue concave bending line
            NSVGpath *path = shape->paths;
            Stroke stroke = addBendingStrokeFromSVG(path);
            concave_strokes.push_back(stroke);
        } else {
            std::cout << "Stroke type not recognized" << std::endl;
        }
    }
    nsvgDelete(image);
}

Sketch::Stroke Sketch::addBendingStrokeFromSVG(NSVGpath *path) {
    Stroke stroke;
    float curr_length = 0;
    bool first_point_of_segment = true;
    bool add_segment_to_stroke = false;
    float max_segment_length = getBendingStrokeSegmentLength();
    Eigen::Vector2f last_point(0,0);
    CurvatureStrokeSegment segment;
    for (int i = 0; i < path->npts - 1; i += 2) {

        // make sure this is not a duplicate of the previous point
        if (i > 0) {
            if (path->pts[i-2] == path->pts[i] && path->pts[i-1] == path->pts[i+1]) {
                std::cout << "found duplicate" << std::endl;
                continue;
            }
        }

        std::shared_ptr<StrokePoint> point = std::make_shared<StrokePoint>();
        point->coordinates = Eigen::Vector2f(path->pts[i], path->pts[i+1]);
        if (!first_point_of_segment) {
            curr_length += (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
        }
        first_point_of_segment = false;
        if (curr_length > max_segment_length) {
            first_point_of_segment = true;
            curr_length = 0;
            add_segment_to_stroke = true;
        }

        Eigen::Vector2f backward_dir(0,0);
        float backward_dist = 1;
        Eigen::Vector2f forward_dir(0,0);
        float forward_dist = 1;
        if (i == 0) {
            forward_dir = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            forward_dist = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
        } else if (i == path->npts - 2) {
            backward_dir = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            backward_dist = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
        } else {
            forward_dir = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            forward_dist = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
            backward_dir = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            backward_dist = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
        }

        point->tangent_dir = (backward_dir / backward_dist + forward_dir / forward_dist).normalized();

        // if this statement executes, then the next point is a duplicate of the current point.
        // We must look ahead two points to compute the tangent directions
        if ((std::isnan(point->tangent_dir.x()) || std::isnan(point->tangent_dir.y())) && i < path->npts - 5 &&
            (path->pts[i+4] != path->pts[i] && path->pts[i+5] != path->pts[i+1])) {
            forward_dir = (Eigen::Vector2f(path->pts[i+4], path->pts[i+5]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            forward_dist = (Eigen::Vector2f(path->pts[i+4], path->pts[i+5]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
            backward_dir = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            backward_dist = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();
            point->tangent_dir = (backward_dir / backward_dist + forward_dir / forward_dist).normalized();
        } else if (std::isnan(point->tangent_dir.x()) || std::isnan(point->tangent_dir.y())) {
            point->tangent_dir = (backward_dir / backward_dist).normalized();
        }

        if (std::isnan(point->tangent_dir.x()) || std::isnan(point->tangent_dir.y())) {
            std::cout << "stop" << std::endl;
        }

        assert(!std::isnan(point->tangent_dir.x()) && !std::isnan(point->tangent_dir.y()));
        assert(!std::isinf(point->tangent_dir.x()) && !std::isinf(point->tangent_dir.y()));

        segment.seg.push_back(point);

        if (add_segment_to_stroke || (i == path->npts - 2)) {
            stroke.segments.push_back(segment);
            add_segment_to_stroke = false;
            segment = CurvatureStrokeSegment();
        }
    }
    return stroke;
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

    for (int line_idx = 0; line_idx < data.convex.size(); line_idx++) {
        Stroke segment = addBendingStrokeFromScribble(data.convex[line_idx]);
        convex_strokes.push_back(std::move(segment));
    }

    for (int line_idx = 0; line_idx < data.concave.size(); line_idx++) {
        Stroke segment = addBendingStrokeFromScribble(data.concave[line_idx]);
        concave_strokes.push_back(std::move(segment));
    }
}

Sketch::Stroke Sketch::addBendingStrokeFromScribble(const std::vector<Eigen::Vector2f> &s) {
    Stroke stroke;
    float curr_length = 0;
    bool first_point_of_segment = true;
    bool add_segment_to_stroke = false;
    float max_segment_length = getBendingStrokeSegmentLength();
    Eigen::Vector2f last_point(0,0);
    CurvatureStrokeSegment segment;
    for (int i = 0; i < s.size(); i++) {

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
        if (!first_point_of_segment) {
            curr_length += (s[i] - s[i-1]).norm();
        }
        first_point_of_segment = false;
        if (curr_length > max_segment_length) {
            first_point_of_segment = true;
            curr_length = 0;
            add_segment_to_stroke = true;
        }

        Eigen::Vector2f backward_dir(0,0);
        float backward_dist = 1;
        Eigen::Vector2f forward_dir(0,0);
        float forward_dist = 1;
        if (i == 0) {
            forward_dir = (s[i+1] - s[i]).normalized();
            forward_dist = (s[i+1] - s[i]).norm();
        } else if (i == s.size() - 1) {
            backward_dir = (s[i] - s[i-1]).normalized();
            backward_dist = (s[i] - s[i-1]).norm();
        } else {
            forward_dir = (s[i+1] - s[i]).normalized();
            forward_dist = (s[i+1] - s[i]).norm();
            backward_dir = (s[i] - s[i-1]).normalized();
            backward_dist = (s[i] - s[i-1]).norm();
        }

        if (backward_dist == 0 && forward_dist > 0) { point->tangent_dir = (forward_dir / forward_dist).normalized(); }
        else if (forward_dist == 0 && backward_dist > 0) { point->tangent_dir = (backward_dir / backward_dist).normalized(); }
        else { point->tangent_dir = (backward_dir / backward_dist + forward_dir / forward_dist).normalized(); }

        // if this statement executes, then the next point is a duplicate of the current point.
        // We must look ahead two points to compute the tangent directions
        if ((std::isnan(point->tangent_dir.x()) || std::isnan(point->tangent_dir.y())) && i < s.size() - 2 && (s[i] - s[i+2]).norm() > 0) {
            forward_dir = (s[i+2] - s[i]).normalized();
            forward_dist = (s[i+2] - s[i]).norm();
            backward_dir = (s[i+2] - s[i]).normalized();
            backward_dist = (s[i+2] - s[i]).norm();
            point->tangent_dir = (backward_dir / backward_dist + forward_dir / forward_dist).normalized();
        } else if (std::isnan(point->tangent_dir.x()) || std::isnan(point->tangent_dir.y())) {
            point->tangent_dir = (backward_dir / backward_dist).normalized();
        }

        // TODO: need to resolve this by smoothing the stroke. Sometimes, we get something like
        // s[0] = (0,0), s[1] = (0,1), s[2] = (0,0), which ends up failing this assert
        assert(point->tangent_dir.norm() > 0);

        segment.seg.push_back(point);
        if (add_segment_to_stroke || (i == s.size() - 1)) {
            //std::cout << "size: " << segment.seg.size() << std::endl;
            //assert(segment.seg.size() > 2);
            stroke.segments.push_back(segment);
            add_segment_to_stroke = false;
            segment = CurvatureStrokeSegment();
        }
    }
    return stroke;
}

void Sketch::mapIntersectedFacesToStrokes(Mesh &mesh) {
    mapIntersectedFacesToStrokesHelper(mesh, concave_strokes);
    mapIntersectedFacesToStrokesHelper(mesh, convex_strokes);
}

void Sketch::mapIntersectedFacesToStrokesHelper(Mesh &mesh, std::vector<Stroke> &strokes) {
    for (int stroke_idx = 0; stroke_idx < strokes.size(); stroke_idx++) {
        int num_points = getLengthOfStrokeInPoints(strokes[stroke_idx]);
        for (int i = 0; i < num_points; i++) {
            std::shared_ptr<StrokePoint> point = getStrokePointByFlattenedIndex(strokes[stroke_idx], i);
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
                            std::shared_ptr<StrokePoint> prev_point = getStrokePointByFlattenedIndex(strokes[stroke_idx], i-1);
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
        }
    }
}
