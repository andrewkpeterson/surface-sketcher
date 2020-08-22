#ifndef OBJWRITER_H
#define OBJWRITER_H

#include "src/pipeline/Mesh.h"
#include <QFile>
#include <QTextStream>
#include "src/pipeline/Sketch.h"

class OBJWriter
{
public:

    static void writeConstraints(Sketch &sketch) {
        QFile discontinuities("constraints.txt");
        if (discontinuities.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&discontinuities);
            char buf[512];

            for (int i = 0; i < sketch.getConvexStrokes().size(); i++) {
                int size = sketch.getLengthOfStrokeInPoints(sketch.getConvexStrokes()[i]);
                for (int j = 0; j < size; j++) {
                    auto point = sketch.getStrokePointByFlattenedIndex(sketch.getConvexStrokes()[i], j);
                    std::sprintf(buf, "%f %f %f %f", point->coordinates.x(), point->coordinates.y(), point->tangent_dir.x(), point->tangent_dir.y());
                    stream << buf << endl;
                }
            }

            for (int i = 0; i < sketch.getConcaveStrokes().size(); i++) {
                int size = sketch.getLengthOfStrokeInPoints(sketch.getConcaveStrokes()[i]);
                for (int j = 0; j < size; j++) {
                    auto point = sketch.getStrokePointByFlattenedIndex(sketch.getConcaveStrokes()[i], j);
                    std::sprintf(buf, "%f %f %f %f", point->coordinates.x(), point->coordinates.y(), point->tangent_dir.x(), point->tangent_dir.y());
                    stream << buf << endl;
                }
            }
        }
    }

    static void writeOBJ(Mesh &mesh, std::string obj_file, std::string direction_field_file) {
        QFile file(obj_file.c_str());
        if (file.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&file);
            char buf[512];
            auto func1 = [&] (std::shared_ptr<Vertex> v) {
                std::sprintf(buf, "v %f %f %f", v->coords.x(), v->coords.y(), v->height);
                stream << buf << endl;
            };

            auto func2 = [&] (std::shared_ptr<Face> f) {
                std::sprintf(buf, "f %d %d %d", f->vertices[0]->index + 1, f->vertices[1]->index + 1, f->vertices[2]->index + 1);
                stream << buf << endl;
            };

            mesh.forEachVertex(func1);
            mesh.forEachTriangle(func2);
        }

        QFile vectors(direction_field_file.c_str());
        if (vectors.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&vectors);
            char buf[512];
            auto func3 = [&] (std::shared_ptr<Face> f) {
                if (f->valid) {
                    std::sprintf(buf, "%f %f %f %f %f %f", f->centroid.x(), f->centroid.y(), f->u.x(), f->u.y(), f->v.x(), f->v.y());
                    stream << buf << endl;
                }
            };

            mesh.forEachTriangle(func3);
        }

        QFile discontinuities("discontinuities.txt");
        if (discontinuities.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&discontinuities);
            char buf[512];
            auto func4 = [&] (std::shared_ptr<Face> f) {
                if (f->discontinuity) {
                    std::sprintf(buf, "%f %f", f->centroid.x(), f->centroid.y());
                    stream << buf << endl;
                }
            };

            mesh.forEachTriangle(func4);
        }

        /*
        QFile z1z2("z1z2.txt");
        if (z1z2.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&z1z2);
            char buf[512];
            auto func4 = [&] (std::shared_ptr<Face> f) {
                std::sprintf(buf, "%f %f %f %f %f %f", f->centroid.x(), f->centroid.y(), f->z1.x(), f->z1.y(), f->z2.x(), f->z2.y());
                stream << buf << endl;
            };

            mesh.forEachTriangle(func4);
        }

        QFile ab("ab.txt");
        if (ab.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&ab);
            char buf[512];
            auto func4 = [&] (std::shared_ptr<Face> f) {
                std::sprintf(buf, "%f %f %f %f %f %f", f->centroid.x(), f->centroid.y(), f->a.x(), f->a.y(), f->b.x(), f->b.y());
                stream << buf << endl;
            };

            mesh.forEachTriangle(func4);
        }
        */

        QFile edges("edge_vertices.txt");
        if (edges.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&edges);
            char buf[512];
            auto func4 = [&] (std::shared_ptr<const Vertex> v) {
                std::sprintf(buf, "%f %f", v->coords.x(), v->coords.y());
                stream << buf << endl;
            };

            mesh.forEachBoundaryVertex(func4);
        }
    }
};

#endif // OBJWRITER_H
