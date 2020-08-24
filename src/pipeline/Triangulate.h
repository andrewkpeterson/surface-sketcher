#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <iostream>
#include <CGAL/Polygon_2.h>
#include <QFile>
#include <QTextStream>
#include <map>
#include <tuple>

struct FaceInfo2
{
  int nesting_level;
  bool in_domain() const {
    return nesting_level%2 == 1;
  }
};

class Mesh;
class Sketch;

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
//typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
//typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::No_intersection_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag   >  CDT;
//typedef CGAL::Exact_predicates_tag                               Itag;
//typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                  Criteria;
typedef CDT::Face_handle                                          Face_handle;
typedef CDT::Vertex_handle                                        Vertex_handle;


struct IntersectionResult {
    int A_stroke_idx;
    int B_stroke_idx;
    int A_point_idx;
    int B_point_idx;
    float dist;
};

struct Neighbor {
    int neighbor_idx;
    bool adjacent_to_index_0_of_neighbor;
};

struct SegmentedStroke {
    std::vector<Eigen::Vector2f> points;
    std::vector<Neighbor> last_index_neighbors;
    std::vector<Neighbor> index_0_neighbors;
};

class Triangulate {
public:
    static void triangulate(Mesh &mesh, Sketch &sketch);

private:
    static void mark_domains(CDT& cdt, std::map<Face_handle, FaceInfo2> &m);
    static void mark_domains(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border, std::map<Face_handle, FaceInfo2> &m);
    static std::vector<IntersectionResult> findIntersections(const std::vector<Eigen::Vector2f> &A, int A_stroke_idx,
                                                             const std::vector<Eigen::Vector2f> &B, int B_stroke_idx);
    static std::vector<IntersectionResult> intersectionExists(const std::vector<Eigen::Vector2f> &A, const std::vector<Eigen::Vector2f> &B);
    static void findPath(const std::vector<SegmentedStroke> &strokes, int first_stroke, std::vector<Eigen::Vector2f> &boundary, std::vector<int> &strokes_sequence);
    static std::vector<Eigen::Vector2f> findPathHelper(const std::vector<SegmentedStroke> &strokes,
                                                       int first_stroke, int current_stroke, std::vector<Eigen::Vector2f> boundary, std::vector<int> &strokes_sequence,
                                                       std::set<int> used_strokes, bool arrived_at_current_stroke_at_index_0, int depth);
    static void addToBoundary(const SegmentedStroke &stroke, std::vector<Eigen::Vector2f> &boundary, bool met_at_idx_0_neighbor);
    static std::vector<SegmentedStroke> segmentBoundaryStrokes(const std::vector<std::vector<Eigen::Vector2f>> &strokes);
    static constexpr float FINENESS = .1; //.1 normally, .2 for coarse
    static constexpr float CRITERIA_PARAM = .001; //.001 normally
    static constexpr float SAME_POINT_THRESHOLD = .001;
    static constexpr float SAME_INTERSECTION_THRESHOLD = .1;
    static constexpr float LINE_TAIL_LENGTH = .2; // this will need to be tuned
};

#endif // TRIANGULATE_H
