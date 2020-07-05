#include "Sketch.h"
#define NANOSVG_IMPLEMENTATION
#define NANOSVG_CPLUSPLUS
#include "lib/nanosvg/src/nanosvg.h"
#include "Mesh.h"
#include <iostream>

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
    Eigen::Vector2f v2(f->vertices[0]->coords);
    Eigen::Vector2f v3(f->vertices[0]->coords);
    A << v1.x(), v2.x(), v3.x(),
         v1.y(), v2.y(), v3.y(),
         1, 1, 1;
    return A.inverse() * x;
}

Sketch::Sketch(std::string svg_file) {
    NSVGimage *image;
    image = nsvgParseFromFile(svg_file.c_str(), "px", 1);

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
            Stroke stroke = addBendingStroke(path);
            convex_strokes.push_back(stroke);
        } else if (shape->stroke.color == 0xFFFF0000) { // blue concave bending line
            NSVGpath *path = shape->paths;
            Stroke stroke = addBendingStroke(path);
            concave_strokes.push_back(stroke);
        } else {
            std::cout << "Stroke type not recognized" << std::endl;
        }
    }
    nsvgDelete(image);
}

Sketch::Stroke Sketch::addBendingStroke(NSVGpath *path) {
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
            forward_dist = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();// + SMOOTH_DIRECTION;
        } else if (i == path->npts - 2) {
            backward_dir = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            backward_dist = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();//+ SMOOTH_DIRECTION;
        } else {
            forward_dir = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            forward_dist = (Eigen::Vector2f(path->pts[i+2], path->pts[i+3]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm(); //+ SMOOTH_DIRECTION;
            backward_dir = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).normalized();
            backward_dist = (Eigen::Vector2f(path->pts[i-2], path->pts[i-1]) - Eigen::Vector2f(path->pts[i], path->pts[i+1])).norm();//+ SMOOTH_DIRECTION;
        }

        point->tangent_dir = (backward_dir / backward_dist + forward_dir / forward_dist).normalized();
        if (std::isnan(point->tangent_dir.x())) {
            std::cout << "stop" << std::endl;
        }
        segment.seg.push_back(point);

        if (add_segment_to_stroke || (i == path->npts - 2)) {
            stroke.push_back(segment);
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
        for (int segment_idx = 0; segment_idx < strokes[stroke_idx].size(); segment_idx++) {
            for (int point_idx = 0; point_idx < strokes[stroke_idx][segment_idx].seg.size(); point_idx++) {
                std::shared_ptr point = strokes[stroke_idx][segment_idx].seg[point_idx];
                auto func = [&] (std::shared_ptr<Face> f) {
                    // do a distance check so that time is not wasted on far-away triangles
                    if ((f->circumcenter - point->coordinates).norm() < getTriangleAndStrokePointDistanceCheck()) {
                        // check if the StrokePoint is inside the triangle, or if one of the two neighboring
                        // line segments crosses the triangle.

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
                            //f->valid = false;
                        } else {
                            if (segment_idx > 0 || point_idx > 0) {
                            // check if previous line segment crossed triangle. There is never any need to check the next line segment ahead of this point.
                                std::shared_ptr<StrokePoint> prev_point;
                                if (point_idx > 0) {
                                    prev_point = strokes[stroke_idx][segment_idx].seg[point_idx-1];
                                } else if (point_idx == 0) {
                                    prev_point = strokes[stroke_idx][segment_idx-1].seg[(strokes[stroke_idx][segment_idx-1].seg.size()-1)];
                                }
                                if (lineInTriangle(point->coordinates, prev_point->coordinates, v1, v2, v3)) {
                                    LineIntersection intersection;
                                    intersection.points = std::make_pair(point, prev_point);
                                    intersection.dir = (point->coordinates - prev_point->coordinates).normalized();
                                    if (face2lines.find(f) != face2lines.end()) {
                                        face2lines[f].push_back(intersection);
                                    } else {
                                        face2lines[f] = std::vector{intersection};
                                    }
                                    //f->valid = false;
                                }
                            }
                        }
                    }
                };

                mesh.forEachTriangle(func);
            }
        }
    }
}
