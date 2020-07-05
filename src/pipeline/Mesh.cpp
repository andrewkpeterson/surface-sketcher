#include "Mesh.h"

float Mesh::calcTriangleArea(Face_handle f) {
    float a = std::sqrt((f->vertex(0)->point() - f->vertex(1)->point()).squared_length());
    float b = std::sqrt((f->vertex(1)->point() - f->vertex(2)->point()).squared_length());
    float c = std::sqrt((f->vertex(0)->point() - f->vertex(2)->point()).squared_length());
    float s = (a + b + c) / 2.0f;
    return std::sqrt(s * (s - a) * (s - b) * (s - c));
}

float Mesh::calcTriangleArea(const Eigen::Vector2f v1, const Eigen::Vector2f v2, const Eigen::Vector2f v3) {
    float a = (v1 - v2).norm();
    float b = (v2 - v3).norm();
    float c = (v1 - v3).norm();
    float s = (a + b + c) / 2.0f;
    return std::sqrt(s * (s - a) * (s - b) * (s - c));
}

// calculate the combined area of the two triangles created by the shared edge and centroids of two adjacent triangles
float Mesh::calcEFGArea(const Face *f, const Face *g) {
    Eigen::Vector2f f_centroid = (f->vertices[0]->coords + f->vertices[1]->coords + f->vertices[2]->coords) / 3.0f;
    Eigen::Vector2f g_centroid = (g->vertices[0]->coords + g->vertices[1]->coords + g->vertices[2]->coords) / 3.0f;

    std::shared_ptr<Vertex> v1 = nullptr;
    std::shared_ptr<Vertex> v2 = nullptr;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (f->vertices[i] == g->vertices[j]) {
                if (v1 == nullptr) {
                    v1 = f->vertices[i];
                } else {
                    v2 = f->vertices[i];
                }
            }
        }
    }

    return calcTriangleArea(v1->coords, v2->coords, f_centroid) + calcTriangleArea(v1->coords, v2->coords, g_centroid);
}

void Mesh::init(std::map<Face_handle, bool> &info) {
    int next_face_index = 0;
    int next_vertex_index = 0;
    int num_faces = 0;
    int num_verts = 0;
    float total_area = 0;

    std::unordered_map<Face_handle, std::shared_ptr<Face>> visited_faces;
    std::unordered_map<Vertex_handle, std::shared_ptr<Vertex>> visited_vertices;

    // get all of the faces and vertices in the mesh
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        if (!cdt.is_infinite(it) && info.find(it)->second) {
            std::shared_ptr<Face> f = std::make_shared<Face>();
            index2face[next_face_index] = f;
            f->index = next_face_index;
            f->circumcenter = calculateCircumcenter(it);
            f->area = calcTriangleArea(it);
            total_area += f->area;
            visited_faces[it] = f;
            next_face_index++;
            num_faces++;

            for (int i = 0; i < 3; i++) {
                if (visited_vertices.find(it->vertex(i)) == visited_vertices.end()) {
                    std::shared_ptr<Vertex> v = std::make_shared<Vertex>();
                    v->coords = Eigen::Vector2f(it->vertex(i)->point().x(), it->vertex(i)->point().y());
                    v->index = next_vertex_index;
                    v->height = 0;
                    visited_vertices[it->vertex(i)] = v;
                    next_vertex_index++;
                    num_verts++;
                }
            }
        }
    }

    // set all of the neighbors and vertices of the faces in the mesh
    for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++) {
        if (!cdt.is_infinite(it) && info.find(it)->second) {
            auto f = visited_faces[it];
            for (int i = 0; i < 3; i++) {
                if (!cdt.is_infinite(it->neighbor(i)) && info.find(it->neighbor(i))->second) {
                    f->neighbors.push_back(visited_faces[it->neighbor(i)].get());
                }
                f->vertices.push_back(visited_vertices[it->vertex(i)]);
            }
        }
    }

    m_num_triangles = num_faces;
    m_total_area = total_area;
    m_num_vertices = num_verts;
}

void Mesh::forEachTriangle(const std::function<void(std::shared_ptr<Face>)> &func) {
    for (auto it = index2face.begin(); it != index2face.end(); it++) {
        func(it->second);
    }
}

void Mesh::forEachConstTriangle(const std::function<void(std::shared_ptr<const Face>)> &func) const {
    for (auto it = index2face.begin(); it != index2face.end(); it++) {
        func(it->second);
    }
}

void Mesh::forEachPairOfNeighboringTriangles(const std::function<void(Face*, Face*)> &func) {
    std::unordered_set<Face*> visited_faces;
    for (auto it = index2face.begin(); it != index2face.end(); it++) {
        for (int i = 0; i < it->second->neighbors.size(); i++) {
            if (visited_faces.find(it->second->neighbors[i]) == visited_faces.end()) {
                func(it->second.get(), it->second->neighbors[i]);
            }
        }
        visited_faces.insert(it->second.get());
    }
}

void Mesh::forEachConstPairOfNeighboringTriangles(const std::function<void(const Face*, const Face*)> &func) const {
    std::unordered_set<Face*> visited_faces;
    for (auto it = index2face.begin(); it != index2face.end(); it++) {
        for (int i = 0; i < it->second->neighbors.size(); i++) {
            if (visited_faces.find(it->second->neighbors[i]) == visited_faces.end()) {
                func(it->second.get(), it->second->neighbors[i]);
            }
        }
        visited_faces.insert(it->second.get());
    }
}

void Mesh::forEachVertex(const std::function<void(std::shared_ptr<Vertex>)> &func) {
    for (auto it = index2vertex.begin(); it != index2vertex.end(); it++) {
        func(it->second);
    }
}

void Mesh::forEachConstVertex(const std::function<void(std::shared_ptr<const Vertex>)> &func) const {
    for (auto it = index2vertex.begin(); it != index2vertex.end(); it++) {
        func(it->second);
    }
}

Eigen::Vector2f Mesh::calculateCircumcenter(const Face_handle f) {
    Point p = CGAL::circumcenter(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
    return Eigen::Vector2f(p.x(), p.y());
}
