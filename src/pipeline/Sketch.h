#ifndef SKETCH_H
#define SKETCH_H

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "Triangulate.h"

struct NSVGpath;

struct Face;

class Sketch
{
public:
    Sketch(std::string svg_file);


    struct StrokePoint {
        Eigen::Vector2f coordinates;
        Eigen::Vector2f tangent_dir;
        std::shared_ptr<Face> triangle = nullptr; // the triangle this point is inside of
        Eigen::Vector3f barycentric_coordinates; // barycentric coordinates of the point inside the triangle
    };

    struct CurvatureStrokeSegment {
        // a segment begins at the terminal point of the previous segment and terminates at its last point.
        // The exception to this is the first segment in a stroke, which begins at its own first point.
        std::vector<std::shared_ptr<StrokePoint>> seg;
        std::set<std::shared_ptr<Face>> faces_intersected;
        float curvature_magnitude;
    };

    using LineIntersection = std::pair<std::shared_ptr<StrokePoint>, std::shared_ptr<StrokePoint>>;

    using Stroke = std::vector<CurvatureStrokeSegment>;
    int width;
    int height;
    int diagonal_of_bounding_box;

    const std::vector<std::vector<Eigen::Vector2f>> &getBoundaryStrokes() const {
        return boundary_strokes;
    }

    Stroke addBendingStroke(NSVGpath *path);
    float getBendingStrokeSegmentLength() { return diagonal_of_bounding_box / 30.0f; }
    float getTriangleAndStrokePointDistanceCheck() const { return diagonal_of_bounding_box / 10.0f; }
    void mapIntersectedFacesToStrokes(Mesh &mesh);

private:
    std::vector<std::vector<Eigen::Vector2f>> boundary_strokes;
    std::vector<Stroke> convex_strokes;
    std::vector<Stroke> concave_strokes;
    // we use shared_ptrs because the StrokePoints have to belong to the strokes and the maps
    std::map<std::shared_ptr<Face>, std::vector<std::shared_ptr<StrokePoint>>> face2strokepoints;
    std::map<std::shared_ptr<Face>, std::vector<LineIntersection>> face2lines;

    void mapIntersectedFacesToStrokesHelper(Mesh &mesh, std::vector<Stroke> &strokes);
};

#endif // SKETCH_H
