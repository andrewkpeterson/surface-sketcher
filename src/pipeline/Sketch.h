#ifndef SKETCH_H
#define SKETCH_H

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "ScribbleArea.h"
#include "Triangulate.h"

struct NSVGpath;

struct Face;

class Sketch
{
public:
    Sketch(std::string svg_file);
    Sketch(const SketchData &data);


    struct StrokePoint {
        Eigen::Vector2f coordinates;
        Eigen::Vector2f tangent_dir;
        std::shared_ptr<Face> triangle = nullptr; // the triangle this point is inside of
        Eigen::Vector3f barycentric_coordinates; // barycentric coordinates of the point inside the triangle
        Eigen::Vector3f coords3d();
        float curvature_value;
        int strokeId;
    };

    struct CurvatureStrokeSegment {
        // a segment begins at the terminal point of the previous segment and terminates at its last point.
        // The exception to this is the first segment in a stroke, which begins at its own first point.
        std::vector<std::shared_ptr<StrokePoint>> seg;
    };

    struct LineIntersection {
        std::pair<std::shared_ptr<StrokePoint>, std::shared_ptr<StrokePoint>> points;
        Eigen::Vector2f dir;
    };

    struct Stroke {
        std::vector<CurvatureStrokeSegment> segments;
        int id;
    };

    int width;
    int height;
    int diagonal_of_bounding_box;

    const std::vector<std::vector<Eigen::Vector2f>> &getBoundaryStrokes() const {
        return boundary_strokes;
    }

    Stroke addBendingStrokeFromSVG(NSVGpath *path);
    Stroke addBendingStrokeFromScribble(const std::vector<Eigen::Vector2f> &s, int id);
    float getBendingStrokeSegmentLength() { return diagonal_of_bounding_box / 30.0f; }
    float getTriangleAndStrokePointDistanceCheck() const { return diagonal_of_bounding_box / 10.0f; }
    void mapIntersectedFacesToStrokes(Mesh &mesh);
    bool checkStrokePoints(std::shared_ptr<Face> f) const { return face2strokepoints.find(f) != face2strokepoints.end(); }
    bool checkStrokeLines(std::shared_ptr<Face> f) const { return face2lines.find(f) != face2lines.end(); }
    const std::vector<std::shared_ptr<StrokePoint>>& getStrokePoints(std::shared_ptr<Face> f) const {
        assert(face2strokepoints.find(f) != face2strokepoints.end());
        return face2strokepoints.at(f);
    }
    const std::vector<LineIntersection>& getStrokeLines(std::shared_ptr<Face> f) const {
        assert(face2lines.find(f) != face2lines.end());
        return face2lines.at(f);
    }
    std::vector<Stroke>& getConvexStrokes() { return convex_strokes; }
    std::vector<Stroke>& getConcaveStrokes() { return concave_strokes; }
    static int getLengthOfStrokeInPoints(const Stroke &stroke);
    static std::shared_ptr<StrokePoint> getStrokePointByFlattenedIndex(const Stroke &stroke, int idx);

    static constexpr double SKETCH_SCALE = 100.0; // factor by which the sketch is scaled down

private:
    std::vector<std::vector<Eigen::Vector2f>> boundary_strokes;
    std::vector<Stroke> convex_strokes;
    std::vector<Stroke> concave_strokes;
    // we use shared_ptrs because the StrokePoints have to belong to the strokes and the maps
    std::map<std::shared_ptr<Face>, std::vector<std::shared_ptr<StrokePoint>>> face2strokepoints;
    std::map<std::shared_ptr<Face>, std::vector<LineIntersection>> face2lines;

    void mapIntersectedFacesToStrokesHelper(Mesh &mesh, std::vector<Stroke> &strokes);

    static constexpr int SVG_SIZE_PARAM = 10;
};

#endif // SKETCH_H
