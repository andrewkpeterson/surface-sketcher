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
                last_successful_idx = a_idx;
                final_results.push_back(r);
            }
        }
    }
    return final_results;
}

void Triangulate::addToBoundary(const SegmentedStroke &stroke, std::vector<Eigen::Vector2f> &boundary, bool met_at_idx_0_of_neighbor) {
    if (met_at_idx_0_of_neighbor) {
        boundary.insert(boundary.end(), stroke.points.begin(), stroke.points.end());
    } else {
        boundary.insert(boundary.end(), stroke.points.rbegin(), stroke.points.rend());
    }
}

std::vector<Eigen::Vector2f> Triangulate::makeBoundary(const std::vector<SegmentedStroke> &strokes, const std::set<int> invalid_strokes) {
    std::vector<Eigen::Vector2f> boundary;
    std::set<int> used_strokes;
    int initial_stroke_idx;
    int current_stroke_idx;

    // find the first stroke that is not invalid and set it as the initial stroke
    for (int i = 0; i < strokes.size(); i++) {
        if (invalid_strokes.find(i) == invalid_strokes.end())  {
            current_stroke_idx = i;
            initial_stroke_idx = current_stroke_idx;
            addToBoundary(strokes[current_stroke_idx], boundary, true);
            used_strokes.insert(current_stroke_idx);
            break;
        }
    }

    while (used_strokes.size() < strokes.size() - invalid_strokes.size()) {
        bool success = false;
        bool met_at_index_0_of_neighbor = false;
        int next;
        const auto &current_stroke = strokes[current_stroke_idx];

        // look for the next stroke to visit
        // only visit an index_0_neighbor if we are not at the intial_stroke_idx
        if (current_stroke_idx != initial_stroke_idx) {
            for (int i = 0; i < current_stroke.index_0_neighbors.size(); i++) {
                int idx = current_stroke.index_0_neighbors[i].neighbor_idx;
                if (used_strokes.find(idx) == used_strokes.end() && invalid_strokes.find(idx) == invalid_strokes.end()) {
                    success = true;
                    next = idx;
                    met_at_index_0_of_neighbor = current_stroke.index_0_neighbors[i].adjacent_to_index_0_of_neighbor;
                }
            }
        }
        for (int i = 0; i < current_stroke.last_index_neighbors.size(); i++) {
            int idx = current_stroke.last_index_neighbors[i].neighbor_idx;
            if (used_strokes.find(idx) == used_strokes.end() && invalid_strokes.find(idx) == invalid_strokes.end()) {
                success = true;
                next = idx;
                met_at_index_0_of_neighbor = current_stroke.last_index_neighbors[i].adjacent_to_index_0_of_neighbor;
            }
        }

        // if we found the next stroke and it is not the initial stroke, then we can add it to the boundary
        // otherwise, we start over at a new stroke
        if (success) {
            addToBoundary(strokes[next], boundary, met_at_index_0_of_neighbor);
            used_strokes.insert(next);
            current_stroke_idx = next;
        } else {
            // add the stroke to the boundary because it might indicate a hole in the surface
            addToBoundary(current_stroke, boundary, true);
            used_strokes.insert(current_stroke_idx);
            for (int i = 0; i < strokes.size(); i++) {
                if (used_strokes.find(i) == used_strokes.end() && invalid_strokes.find(i) == invalid_strokes.end())  {
                    current_stroke_idx = i;
                    initial_stroke_idx = i;
                }
            }
        }
    }

    return boundary;
}

std::pair<std::vector<SegmentedStroke>, std::set<int>> Triangulate::segmentBoundaryStrokes(const std::vector<std::vector<Eigen::Vector2f>> &strokes) {
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
                    segmented.points = {strokes[stroke_idx][i]};
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
                if ((A[0] - B[0]).norm() < CONNECTION_THRESHOLD) {
                    segmented_strokes[i].index_0_neighbors.push_back({j, true});
                }
                if ((A[0] - B[B.size() - 1]).norm() < CONNECTION_THRESHOLD) {
                    segmented_strokes[i].index_0_neighbors.push_back({j, false});
                }
                if ((A[A.size() - 1] - B[0]).norm() < CONNECTION_THRESHOLD) {
                    segmented_strokes[i].last_index_neighbors.push_back({j, true});
                }
                if ((A[A.size() - 1] - B[B.size() - 1]).norm() < CONNECTION_THRESHOLD) {
                    segmented_strokes[i].last_index_neighbors.push_back({j, false});
                }
            }
        }
    }

    // mark any segmented strokes with neighbors only on one side
    std::set<int> invalid_strokes;
    for (int i = 0; i < segmented_strokes.size(); i++) {
        int size1 = segmented_strokes[i].last_index_neighbors.size();
        int size2 = segmented_strokes[i].index_0_neighbors.size();
        if ((size1 == 0 && size2 > 0) || (size1 > 0 && size2 == 0)) {
            invalid_strokes.insert(i);
        }
    }

    QFile magnitudes("boundary.txt");
    if (magnitudes.open(QIODevice::ReadWrite | QFile::Truncate)) {
        QTextStream stream(&magnitudes);
        char buf[512];
        for (int i = 0; i < segmented_strokes.size(); i++) {
            //if (invalid_strokes.find(i) == invalid_strokes.end()) {
                for (int j = 0; j < segmented_strokes[i].points.size(); j++) {
                    std::sprintf(buf, "%f %f %d", segmented_strokes[i].points[j].x(), segmented_strokes[i].points[j].y(), i);
                    stream << buf << endl;
                }
            //}
        }
    }

    return std::pair(segmented_strokes, invalid_strokes);
}

void Triangulate::triangulate(Mesh &mesh, Sketch &sketch) {

    const auto &strokes = sketch.getBoundaryStrokes();
    std::vector<std::vector<Eigen::Vector2f>> boundary_strokes;
    for (int i = 0; i < strokes.size(); i++) { boundary_strokes.push_back(strokes[i].points); }
    const auto p = segmentBoundaryStrokes(boundary_strokes);
    std::vector<Eigen::Vector2f> boundary = makeBoundary(p.first, p.second);

    // only add every third boundary point to the polygon so that the constraints do not
    // slow the program down
    Polygon_2 polygon;
    for (int i = 0; i < boundary.size(); i += NUMBER_OF_BOUNDARY_POINTS_TO_SKIP) {
        polygon.push_back(Point(boundary[i].x(), boundary[i].y()));
        std::cout << "x: " << boundary[i].x() << " y: " << boundary[i].y() << std::endl;
    }

    QFile magnitudes("boundary.txt");
    if (magnitudes.open(QIODevice::ReadWrite | QFile::Truncate)) {
        QTextStream stream(&magnitudes);
        char buf[512];
        for (int i = 0; i < boundary.size(); i += NUMBER_OF_BOUNDARY_POINTS_TO_SKIP) {
            std::sprintf(buf, "%f %f", boundary[i].x(), boundary[i].y());
            stream << buf << endl;
        }
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

    PolygonMesh pmesh;
    std::map<Vertex_handle, PolygonMesh::vertex_index> vertices;
    std::set<PolygonMesh::face_index> faces;
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        if (!cdt.is_infinite(it) && domain_info[it].in_domain()) {
            for (int i = 0; i < 3; i++) {
                if (vertices.find(it->vertex(i)) == vertices.end()) {
                    vertices[it->vertex(i)] = pmesh.add_vertex(Kernel::Point_3(it->vertex(i)->point().x(), it->vertex(i)->point().y(), 0));
                }
            }
            faces.insert(pmesh.add_face(vertices[it->vertex(0)], vertices[it->vertex(1)], vertices[it->vertex(2)]));
        }
    }

    //CGAL::Subdivision_method_3::Loop_subdivision(pmesh, CGAL::parameters::number_of_iterations(2));
    std::cout << "done" << std::endl;

    mesh.addSurfaceToMesh(info, sketch);
    mesh.squishTriangulation();

    /*
    // save the 2D triangulation to an obj so that it can be checked
    QFile file("triangulation.obj");
    if (file.open(QIODevice::ReadWrite | QFile::Truncate)) {
        QTextStream stream(&file);
        char buf[512];
        int count = 1;
        std::map<PolygonMesh::Vertex_index, int> m;
        for (PolygonMesh::Vertex_index idx : pmesh.vertices()) {
            std::sprintf(buf, "v %f %f %f", pmesh.point(idx).x(), pmesh.point(idx).y(), pmesh.point(idx).z());
            stream << buf << endl;
            m[idx] = count;
            count++;
        }

        for (PolygonMesh::Face_index idx : pmesh.faces()) {
            stream << "f ";
            CGAL::Vertex_around_face_iterator<PolygonMesh> vbegin, vend;
            for(PolygonMesh::Vertex_index vd : pmesh.vertices_around_face(pmesh.halfedge(idx))) {
                std::sprintf(buf, "%d ", int(m[vd]));
                stream << buf;
            }

            stream << endl;
        }
    }
    */
}
