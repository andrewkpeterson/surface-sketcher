#ifndef OBJWRITER_H
#define OBJWRITER_H

#include "src/pipeline/Mesh.h"
#include <QFile>
#include <QTextStream>

class OBJWriter
{
public:
    static void writeOBJ(Mesh &mesh, std::string obj_file, std::string direction_field_file) {
        QFile file(obj_file.c_str());
        if (file.open(QIODevice::ReadWrite | QFile::Truncate)) {
            QTextStream stream(&file);
            char buf[512];
            auto func1 = [&] (std::shared_ptr<Face> f) {
                if (f->valid) {
                    std::sprintf(buf, "v %f %f %f", f->vertices[0]->coords.x(), f->vertices[0]->coords.y(), f->vertices[0]->height);
                    stream << buf << endl;
                    std::sprintf(buf, "v %f %f %f", f->vertices[1]->coords.x(), f->vertices[1]->coords.y(), f->vertices[1]->height);
                    stream << buf << endl;
                    std::sprintf(buf, "v %f %f %f", f->vertices[2]->coords.x(), f->vertices[2]->coords.y(), f->vertices[2]->height);
                }
            };

            int count = 1;
            auto func2 = [&] (std::shared_ptr<Face> f) {
                if (f->valid) {
                    std::sprintf(buf, "f %d//%d %d//%d %d//%d", count, count, count + 1, count + 1, count + 2, count + 2);
                    stream << buf << endl;
                    count += 3;
                }
            };

            mesh.forEachTriangle(func1);
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
    }
};

#endif // OBJWRITER_H
