#include "Triangulate.h"
#include <CGAL/lloyd_optimize_mesh_2.h>
#include "Mesh.h"
#include "Sketch.h"


void Triangulate::mark_domains(CDT &ct, Face_handle start, int index, std::list<CDT::Edge> &border, std::map<Face_handle, FaceInfo2> &m)
{
  if(m[start].nesting_level != -1){
    return;
  }
  std::list<Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    Face_handle fh = queue.front();
    queue.pop_front();
    if(m[fh].nesting_level == -1){
      m[fh].nesting_level = index;
      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        Face_handle n = fh->neighbor(i);
        if(m[n].nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

void Triangulate::mark_domains(CDT &cdt, std::map<Face_handle, FaceInfo2> &m)
{
  for(CDT::Face_handle f : cdt.all_face_handles()){
    m[f].nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border, m);
  while(!border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    Face_handle n = e.first->neighbor(e.second);
    if(m[n].nesting_level == -1){
      mark_domains(cdt, n, m[e.first].nesting_level+1, border, m);
    }
  }
}

void removeFace(CDT &cdt, std::map<Face_handle, FaceInfo2> &domain_info, Face_handle f);

void checkFace(CDT &cdt, std::map<Face_handle, FaceInfo2> &domain_info, Face_handle f) {
    if (!cdt.is_infinite(f) && domain_info[f].in_domain()) {
        int count = 0;
        if (!cdt.is_infinite(f->neighbor(0)) && domain_info[f->neighbor(0)].in_domain() && Mesh::calcTriangleArea(f->neighbor(0)) > .0001) { count++; }
        if (!cdt.is_infinite(f->neighbor(1)) && domain_info[f->neighbor(1)].in_domain() && Mesh::calcTriangleArea(f->neighbor(1)) > .0001) { count++; }
        if (!cdt.is_infinite(f->neighbor(2)) && domain_info[f->neighbor(2)].in_domain() && Mesh::calcTriangleArea(f->neighbor(2)) > .0001) { count++; }
        // if the number of valid (not infinite, in domain, not too small) neighboring triangles is less than 2, or if the triangle is just really
        // small, then we should get rid of it an consider getting rid of its neighboring triangles as well
        if (count < 2 || Mesh::calcTriangleArea(f) < .0001) {
            removeFace(cdt, domain_info, f);
        }
    }
}

void removeFace(CDT &cdt, std::map<Face_handle, FaceInfo2> &domain_info, Face_handle f) {
    domain_info[f].nesting_level = 0; // this will mark the face as outside of the domain, effectively deleting it
    // check the neighboring triangles if they are not infinite are currently included in the final mesh
    if (!cdt.is_infinite(f->neighbor(0)) && domain_info[f->neighbor(0)].in_domain()) { checkFace(cdt, domain_info, f->neighbor(0)); }
    if (!cdt.is_infinite(f->neighbor(1)) && domain_info[f->neighbor(1)].in_domain()) { checkFace(cdt, domain_info, f->neighbor(1)); }
    if (!cdt.is_infinite(f->neighbor(2)) && domain_info[f->neighbor(2)].in_domain()) { checkFace(cdt, domain_info, f->neighbor(2)); }
}

void removeArtifacts(CDT &cdt, std::map<Face_handle, FaceInfo2> &domain_info) {
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        checkFace(cdt, domain_info, it);
    }
}


void Triangulate::triangulate(Mesh &mesh, const Sketch &sketch) {
    Polygon_2 polygon;

    // turn the boundary lines into a polygon that we can pass to CGAL
    // TODO: Process lines so that the entire boundary is a single line with no
    // separate intersecting segments. This will result in a better triangulation
    // without any bowtie-like shapes.
    /*
    for (int i = 0; i < sketch.getBoundaryStrokes().size(); i++) {

        // ensure that the boundary is closed by connecting the later endpoint of this line to the
        // closest endpoint of another boundary line
        Eigen::Vector2f curr_endpoint = sketch.getBoundaryStrokes()[i][sketch.getBoundaryStrokes()[i].size() - 1];
        float min_distance = INFINITY;
        Eigen::Vector2f best_endpoint;
        for (int j = 0; j < sketch.getBoundaryStrokes().size(); j++) {
            float distance = (sketch.getBoundaryStrokes()[j][0] - curr_endpoint).norm();
            if (distance > 0 && distance < min_distance) {
                min_distance = distance;
                best_endpoint = sketch.getBoundaryStrokes()[j][0];
            }
            distance = (sketch.getBoundaryStrokes()[j][sketch.getBoundaryStrokes()[j].size() - 1] - curr_endpoint).norm();
            if (distance > 0 && distance < min_distance) {
                min_distance = distance;
                best_endpoint = sketch.getBoundaryStrokes()[j][sketch.getBoundaryStrokes()[j].size() - 1];
            }
        }

        // once the best endpoint is found, make sure that the last 5 points
        // on the line are not closer for some reason

        for (int j = 0; j < sketch.getBoundaryStrokes()[i].size() - 2; j++) {
            polygon.push_back(Point(sketch.getBoundaryStrokes()[i][j].x(), sketch.getBoundaryStrokes()[i][j].y()));
        }


        polygon.push_back(Point(best_endpoint.x(), best_endpoint.y()));
    }
    */

    // test polygon
    polygon.push_back(Point(0,0));
    polygon.push_back(Point(0,100));
    polygon.push_back(Point(100,100));
    polygon.push_back(Point(100,0));


    CDT &cdt = mesh.getCDT();

    cdt.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);

    CGAL::refine_Delaunay_mesh_2(cdt, Criteria(CRITERIA_PARAM, FINENESS));

    // mark which faces are actually inside the constraints. If we don't do this step
    // then the entire convex hull will be triangulated.
    std::map<Face_handle, FaceInfo2> domain_info;
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        domain_info[it] = FaceInfo2();
    }
    mark_domains(cdt, domain_info);

    // NOTE: a face is valid if it is within the constraints and it is not an infinite face

    // get rid of artifacts by recursively getting rid of all faces that have 1 neighbor.
    // this will result in a mesh with much fewer artifacts
    //removeArtifacts(cdt, domain_info); **********************************************************

    // initialize mesh struct
    std::map<Face_handle, bool> info;
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        if (!cdt.is_infinite(it) && domain_info[it].in_domain()) {
            info[it] = true;
        } else {
            info[it] = false;
        }
    }
    mesh.init(info);

    // save the 2D triangulation to an obj so that it can be checked
    QFile file("triangulation.obj");
    if (file.open(QIODevice::ReadWrite | QFile::Truncate)) {
        QTextStream stream(&file);
        char buf[512];
        for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
            if (domain_info[it].in_domain()) {
                std::sprintf(buf, "v %f %f %f", it->vertex(0)->point().x(), it->vertex(0)->point().y(), 0.0f);
                stream << buf << endl;
                std::sprintf(buf, "v %f %f %f", it->vertex(1)->point().x(), it->vertex(1)->point().y(), 0.0f);
                stream << buf << endl;
                std::sprintf(buf, "v %f %f %f", it->vertex(2)->point().x(), it->vertex(2)->point().y(), 0.0f);
                stream << buf << endl;
            }
        }
        int count = 1;
        for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
            if (domain_info[it].in_domain()) {
                std::sprintf(buf, "f %d %d %d", count, count + 1, count + 2);
                stream << buf << endl;
                count += 3;
            }
        }
    }
}
