#ifndef MESH_H
#define MESH_H

#include <vector>
#include <Eigen/Dense>
#include "Triangulate.h"
#include <functional>
#include <tuple>

class Mesh;

struct Vertex {
    int index;
    Eigen::Vector2f coords;
    float height = 0;

    Eigen::Vector3f coords3d() {
        return Eigen::Vector3f(coords[0], coords[1], height);
    }
};

struct Face {
    int index; // the index of this face for matrices and vectors
    bool valid = true;

    Eigen::Vector2f circumcenter;
    float area;

    std::vector<Face*> neighbors; // use raw pointers to neighbors to avoid shared ptr cycles
    std::vector<std::shared_ptr<Vertex>> vertices;

    Eigen::Vector2f u = Eigen::Vector2f(0,0);
    Eigen::Vector2f v = Eigen::Vector2f(0,0);
    float lambda_u = 0;
    float lambda_v = 0;

    Eigen::Vector3f normal() {

        return (vertices[0]->coords3d() - vertices[1]->coords3d()).cross(vertices[0]->coords3d() - vertices[2]->coords3d()).normalized();
    }
};

class Mesh
{   
public:
    void init(std::map<Face_handle, bool> &info);

    void forEachTriangle(const std::function<void(std::shared_ptr<Face>)> &func);
    void forEachPairOfNeighboringTriangles(const std::function<void(Face*, Face*)> &func);
    void forEachVertex(const std::function<void(std::shared_ptr<Vertex>)> &func);

    void forEachConstTriangle(const std::function<void(std::shared_ptr<const Face>)> &func) const;
    void forEachConstPairOfNeighboringTriangles(const std::function<void(const Face*, const Face*)> &func) const;
    void forEachConstVertex(const std::function<void(std::shared_ptr<const Vertex>)> &func) const;

    static float calcEFGArea(const Face *f, const Face *g);
    int getNumTriangles() const { return m_num_triangles; }
    int getNumVertices() const { return m_num_vertices; }
    float getTotalArea() const { return m_total_area; }
    CDT &getCDT() { return cdt; } // this function is necessary to share the cdt with the triangulation class

    std::shared_ptr<const Face> getConstFace(int i) { return index2face[i]; }
    std::shared_ptr<Face> getFace(int i) { return index2face[i]; }

    static float calcTriangleArea(const Face_handle f);
    static float calcTriangleArea(const Eigen::Vector2f v1, const Eigen::Vector2f v2, const Eigen::Vector2f v3);

private:
    Eigen::Vector2f calculateCircumcenter(const Face_handle f);

    CDT cdt; // constrained delaunay triangulation

    std::unordered_map<int, std::shared_ptr<Face>> index2face;
    std::unordered_map<int, std::shared_ptr<Vertex>> index2vertex;

    int m_num_triangles;
    int m_num_vertices;
    float m_total_area;

};

#endif // MESH_H
