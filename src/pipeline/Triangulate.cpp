#include "Triangulate.h"
#include <CGAL/lloyd_optimize_mesh_2.h>
#include "Mesh.h"
#include "Sketch.h"

// Both of these "mark domains" methods are helpers that allow us to create a triangulation defined by the
// boundary curves, rather than the convex hull of the points defining the boundary curves. The code for the "mark domains"
// methods is from https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2polygon_triangulation_8cpp-example.html
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

//from https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2polygon_triangulation_8cpp-example.html, as explained above
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

// remove artifacts in the triangulation along the border by recursively removing any faces with less
// than 2 neighbors and removing any faces with very small area
void removeArtifacts(CDT &cdt, std::map<Face_handle, FaceInfo2> &domain_info) {
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        checkFace(cdt, domain_info, it);
    }
}

std::vector<IntersectionResult> Triangulate::findIntersections(const std::vector<Eigen::Vector2f> &A, int A_stroke_idx,
                                                               const std::vector<Eigen::Vector2f> &B, int B_stroke_idx) {
    std::vector<IntersectionResult> final_results;
    bool found_one = false;
    int last_successful_idx = 0;
    for (int a_idx = 0; a_idx < A.size(); a_idx++) {
        for (int b_idx = 0; b_idx < B.size(); b_idx++) {
            float dist = (A[a_idx] - B[b_idx]).norm();
            if (dist < SAME_POINT_THRESHOLD && (!found_one || (found_one && (A[a_idx] - A[last_successful_idx]).norm() > SAME_INTERSECTION_THRESHOLD))) {
                IntersectionResult r;
                r.A_stroke_idx = A_stroke_idx;
                r.B_stroke_idx = B_stroke_idx;
                r.A_point_idx = a_idx;
                r.B_point_idx = b_idx;
                r.dist = dist;
                found_one = true;
            }
        }
    }
    return final_results;
}

void Triangulate::addToBoundary(const SegmentedStroke &stroke, std::vector<Eigen::Vector2f> &boundary, bool met_at_idx_0_neighbor) {
    if (met_at_idx_0_neighbor) {
        boundary.insert(boundary.begin(), stroke.points.begin(), stroke.points.end());
    } else {
        boundary.insert(boundary.begin(), stroke.points.rbegin(), stroke.points.rend());
    }
}

std::vector<Eigen::Vector2f> Triangulate::makeBoundary(const std::vector<SegmentedStroke> &strokes) {
    std::set<int> used_strokes;
    auto &initial_stroke = strokes[0];
    int steps_away_from_initial = 0;
    while (used_strokes.size() < strokes.size()) {

    }
}

std::vector<SegmentedStroke> Triangulate::segmentBoundaryStrokes(const std::vector<std::vector<Eigen::Vector2f>> &strokes) {
    // go along each stroke and see if it is intersected at each point
    std::vector<IntersectionResult> results;
    for (int i = 0; i < strokes.size(); i++) {
        for (int j = i + 1; j < strokes.size(); j++) {
            auto results_to_add = findIntersections(strokes[i],i, strokes[j], j);
            results.insert(results.end(), results_to_add.begin(), results_to_add.end());
        }
    }

    std::vector<SegmentedStroke> segmented_strokes;
    for (int stroke_idx = 0; stroke_idx < strokes.size(); stroke_idx++) {
        SegmentedStroke segmented;
        for (int i = 0; i < strokes[stroke_idx].size(); i++) {
            segmented.points.push_back(strokes[stroke_idx][i]);
            for (int result_idx = 0; result_idx < results.size(); result_idx++) {
                bool at_intersection_A = (stroke_idx == results[result_idx].A_stroke_idx && i == results[result_idx].A_point_idx);
                bool at_intersection_B = (stroke_idx == results[result_idx].B_stroke_idx && i == results[result_idx].B_point_idx);
                if (at_intersection_A || at_intersection_B) {
                    segmented_strokes.push_back(segmented);
                    segmented = SegmentedStroke();
                }
            }
        }
        if (segmented.points.size() > 0) {
            segmented_strokes.push_back(segmented);
            segmented = SegmentedStroke();
        }
    }

    // go through the segmented strokes and set their neighbors
    for (int i = 0; i < segmented_strokes.size(); i++) {
        for (int j = 0; j < segmented_strokes.size(); j++) {
            if (j != i) {
                const auto& A = segmented_strokes[i].points;
                const auto& B = segmented_strokes[j].points;
                if ((A[0] - B[0]).norm() < SAME_POINT_THRESHOLD) {
                    segmented_strokes[i].index_0_neighbors.push_back({j, true});
                } else if ((A[0] - B[B.size() - 1]).norm() < SAME_POINT_THRESHOLD) {
                    segmented_strokes[i].index_0_neighbors.push_back({j, false});
                } else if ((A[A.size() - 1] - B[0]).norm() < SAME_POINT_THRESHOLD) {
                    segmented_strokes[i].last_index_neighbors.push_back({j, true});
                } else if ((A[A.size() - 1] - B[B.size() - 1]).norm() < SAME_POINT_THRESHOLD) {
                    segmented_strokes[i].last_index_neighbors.push_back({j, false});
                }
            }
        }
    }

    // remove any segmented strokes with neighbors only on one side
    int i = 0;
    while (i < segmented_strokes.size()) {
        if (segmented_strokes[i].last_index_neighbors.size() == 0 || segmented_strokes[i].index_0_neighbors.size() == 0) {
            segmented_strokes.erase(segmented_strokes.begin() + i);
        } else {
            i++;
        }
    }

    return segmented_strokes;
}

void Triangulate::triangulate(Mesh &mesh, Sketch &sketch) {

    const auto &strokes = sketch.getBoundaryStrokes();
    const auto segmented_boundary_strokes = segmentBoundaryStrokes(strokes);
    std::vector<Eigen::Vector2f> boundary = makeBoundary(segmented_boundary_strokes);
    Polygon_2 polygon;
    for (int i = 0; i < boundary.size(); i++) {
        polygon.push_back(Point(boundary[i].x(), boundary[i].y()));
    }

    // get the length of the boundary
    float length = 0;
    for (int i = 0; i < boundary.size() - 1; i++) {
        length += (boundary[i] - boundary[i+1]).norm();
    }
    sketch.setBoundaryLength(length);

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
    removeArtifacts(cdt, domain_info);

    // initialize mesh struct
    std::map<Face_handle, bool> info;
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        if (!cdt.is_infinite(it) && domain_info[it].in_domain()) {
            info[it] = true;
        } else {
            info[it] = false;
        }
    }
    mesh.addSurfaceToMesh(info);
}
