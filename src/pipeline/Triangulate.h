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



class Triangulate {
public:
    static void triangulate(Mesh &mesh, const Sketch &sketch);

private:
    static void mark_domains(CDT& cdt, std::map<Face_handle, FaceInfo2> &m);
    static void mark_domains(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border, std::map<Face_handle, FaceInfo2> &m);
    static constexpr float FINENESS = 10; //.1 normally, 50 for really coarse, 100 for really really coarse
    static constexpr float CRITERIA_PARAM = .001; //.125 normally, .001 for really coarse
};

#endif // TRIANGULATE_H
