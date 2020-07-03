#ifndef DIRECTIONFIELDSOLVER_H
#define DIRECTIONFIELDSOLVER_H

class Mesh;
class Sketch;

class DirectionFieldSolver
{
public:
    static void initializeDirectionField(Mesh &mesh, Sketch &sketch);
    static void optimizeBendFieldEnergy(Mesh &mesh, Sketch &sketch);
};

#endif // DIRECTIONFIELDSOLVER_H
