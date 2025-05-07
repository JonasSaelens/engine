//
// Created by jonas on 5/7/2025.
//

#ifndef _3D_FIGURES_H
#define _3D_FIGURES_H
#include "structs.h"
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class _3D_Figures {
public:
    static Figure createCube();
    static Figure createTetrahedron();
    static Figure createOctahedron();
    static Figure createIcosahedron();
    static Figure createDodecahedron();
    static Figure createSphere(int n);
    static Figure createCone(int n, double h);
    static Figure createCylinder(int n, double h);
    static Figure createTorus(double r,double R,int n,int m);

};



#endif //_3D_FIGURES_H
