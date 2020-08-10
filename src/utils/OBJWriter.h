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
    }
};

#endif // OBJWRITER_H
