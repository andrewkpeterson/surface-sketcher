#ifndef MESH_H
#define MESH_H

#include <vector>
#include <Eigen/Dense>
#include "Triangulate.h"
#include <functional>
#include <tuple>

class Mesh;
struct Face;

struct Vertex {
    int index;
    Eigen::Vector2f coords;
    float height = 0;
    bool boundary = false;

    Eigen::Vector3f coords3d() {
        return Eigen::Vector3f(coords[0], coords[1], height);
    }
};

struct Face {
    int index; // the index of this face for matrices and vectors
    bool valid = true;

    Eigen::Vector2f circumcenter;
    Eigen::Vector2f centroid;
    float area;

    std::vector<Face*> neighbors; // use raw pointers to neighbors to avoid shared ptr cycles
    std::vector<std::shared_ptr<Vertex>> vertices;

    Eigen::Vector2f u = Eigen::Vector2f::Random().normalized(); //Eigen::Vector2f(0,0);
    Eigen::Vector2f v = Eigen::Vector2f::Random().normalized(); //Eigen::Vector2f(0,0);
    float lambda_u = 0;
    float lambda_v = 0;
    bool discontinuity = false;
    Eigen::Vector2f z1 = Eigen::Vector2f(0,0);
    Eigen::Vector2f z2 = Eigen::Vector2f(0,0);
    Eigen::Vector2f a = Eigen::Vector2f(0,0);
    Eigen::Vector2f b = Eigen::Vector2f(0,0);


    Eigen::Vector3f normal() {
        Eigen::Vector3f n = (vertices[0]->coords3d() - vertices[1]->coords3d()).cross(vertices[0]->coords3d() - vertices[2]->coords3d()).normalized();
        n = n.dot(Eigen::Vector3f(0,0,1)) > 0 ? n : -n;
        return n;
    }
};

class Mesh
{
public:
    void addSurfaceToMesh(std::map<Face_handle, bool> &info);

    void forEachTriangle(const std::function<void(std::shared_ptr<Face>)> &func);
    void forEachPairOfNeighboringTriangles(const std::function<void(Face*, Face*)> &func);
    void forEachVertex(const std::function<void(std::shared_ptr<Vertex>)> &func);

    void forEachTriangle(const std::function<void(std::shared_ptr<const Face>)> &func) const;
    void forEachPairOfNeighboringTriangles(const std::function<void(const Face*, const Face*)> &func) const;
    void forEachVertex(const std::function<void(std::shared_ptr<const Vertex>)> &func) const;
    void forEachBoundaryVertex(const std::function<void(std::shared_ptr<const Vertex>)> &func) const;

    static float calcEFGArea(const Face *f, const Face *g);
    int getNumTriangles() const { return m_num_triangles; }
    int getNumVertices() const { return m_num_vertices; }
    float getTotalArea() const { return m_total_area; }
    CDT &getCDT() { return cdt; } // this function is necessary to share the cdt with the triangulation class

    std::shared_ptr<const Face> getConstFace(int i) { return index2face[i]; }
    std::shared_ptr<Face> getFace(int i) { return index2face[i]; }
    const std::set<Face*> &getEdgeFaces() { return edge_faces; }
    const std::set<std::shared_ptr<Vertex>> &getEdgeVertices() const { return edge_vertices; }

    static float calcTriangleArea(const Face_handle f);
    static float calcTriangleArea(const Eigen::Vector2f v1, const Eigen::Vector2f v2, const Eigen::Vector2f v3);

private:
    Eigen::Vector2f calculateCircumcenter(const Face_handle f);

    CDT cdt; // constrained delaunay triangulation

    std::map<int, std::shared_ptr<Face>> index2face;
    std::map<int, std::shared_ptr<Vertex>> index2vertex;
    std::set<Face*> edge_faces;
    std::set<std::shared_ptr<Vertex>> edge_vertices;

    int m_num_triangles = 0;
    int m_num_vertices = 0;
    float m_total_area = 0;
    int next_face_index = 0;
    int next_vertex_index = 0;

    float m_boundary_length;

};

#endif // MESH_H
