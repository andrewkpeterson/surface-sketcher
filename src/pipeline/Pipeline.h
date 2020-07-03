#ifndef PIPELINE_H
#define PIPELINE_H

#include <string>
#include "Triangulate.h"
#include "Sketch.h"
#include "Mesh.h"
#include "src/utils/OBJWriter.h"
#include "DirectionFieldSolver.h"
#include <iostream>

class Pipeline
{
public:

    static void runPipeline(std::string svg_file, std::string obj_file) {

        std::cout << "Reading in sketch..." << std::endl;
        Sketch sketch(svg_file);
        std::cout << "DONE" << std::endl;

        std::cout << "Triangulating domain..." << std::endl;
        Mesh mesh;
        Triangulate::triangulate(mesh, sketch);
        std::cout << "DONE" << std::endl;

        std::cout << "Processing triangulation and bending lines..." << std::endl;
        sketch.mapIntersectedFacesToStrokes(mesh);
        std::cout << "DONE" << std::endl;

        std::cout << "Initializing curvature direction field..." << std::endl;
        DirectionFieldSolver::initializeDirectionField(mesh, sketch);
        std::cout << "DONE" << std::endl;



        std::cout << "Writing result mesh to file" << std::endl;
        OBJWriter::writeOBJ(mesh, obj_file);
    }

};

#endif // PIPELINE_H
