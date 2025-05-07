//
// Created by jonas on 5/7/2025.
//

#ifndef _3D_MATRIXFUNCTIONS_H
#define _3D_MATRIXFUNCTIONS_H
#include "structs.h"
#include "vector3d.h"


class _3D_MatrixFunctions {
public:
    static Matrix Scale(double scale);
    static Matrix rotateX(double angle);
    static Matrix rotateY(double angle);

    static Matrix rotateZ(double angle);

    static Matrix translate(const Vector3D &vector);
    static void applyTransformation(Figure &fig, const Matrix &m);
    static void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
    static Matrix eyePointTrans(const Vector3D &eyepoint);
    static Point2D doProjectionPoint(const Vector3D &point, double d);
    static Lines2D doProjectionLines(const Figures3D &figures);
};



#endif //_3D_MATRIXFUNCTIONS_H
