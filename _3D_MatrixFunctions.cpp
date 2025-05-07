//
// Created by jonas on 5/7/2025.
//

#include "_3D_MatrixFunctions.h"
#include "_3D_Figures.h"


Matrix _3D_MatrixFunctions::Scale(const double scale) {
        Matrix matrix;
        matrix(1,1) = scale;
        matrix(2,2) = scale;
        matrix(3,3) = scale;
        return matrix;
}
Matrix _3D_MatrixFunctions::rotateX(const double angle) {
        Matrix matrix;
        double radiant = angle * (M_PI /180);
        matrix(2,2) = cos(radiant);
        matrix(2,3) = sin(radiant);
        matrix(3,2) = -sin(radiant);
        matrix(3,3) = cos(radiant);
        return matrix;
}
Matrix _3D_MatrixFunctions::rotateY(const double angle) {
        Matrix matrix;
        double radiant = angle * (M_PI /180);
        matrix(1,1) = cos(radiant);
        matrix(1,3) = -sin(radiant);
        matrix(3,1) = sin(radiant);
        matrix(3,3) = cos(radiant);
        return matrix;
}
Matrix _3D_MatrixFunctions::rotateZ(const double angle) {
        Matrix matrix;
        double radiant = angle * (M_PI /180);
        matrix(1,1) = cos(radiant);
        matrix(1,2) = sin(radiant);
        matrix(2,1) = -sin(radiant);
        matrix(2,2) = cos(radiant);
        return matrix;
}
Matrix _3D_MatrixFunctions::translate(const Vector3D &vector) {
        Matrix matrix;
        matrix(4,1) = vector.x;
        matrix(4,2) = vector.y;
        matrix(4,3) = vector.z;
        return matrix;
}
void _3D_MatrixFunctions::applyTransformation(Figure &fig, const Matrix &m) {
        for (auto &i : fig.points) {
                i *= m;
        }
}
void _3D_MatrixFunctions::toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
        r = sqrt(pow(point.x,2)+pow(point.y,2)+pow(point.z,2));
        theta = atan2(point.y,point.x);
        phi = acos(point.z/r);
}
Matrix _3D_MatrixFunctions::eyePointTrans(const Vector3D &eyepoint) {
        double theta, phi, r;
        toPolar(eyepoint, theta, phi, r);

        Matrix V;

        V(1,1)= -sin(theta);
        V(1,2)= -cos(theta) * cos(phi);
        V(1,3)= cos(theta) * sin(phi);
        V(2,1)= cos(theta);
        V(2,2)= -sin(theta) * cos(phi);
        V(2,3)= sin(theta) * sin(phi);
        V(3,2)= sin(phi);
        V(3,3)= cos(phi);
        V(4,3)= -r;

        return V;
}
Point2D _3D_MatrixFunctions::doProjectionPoint(const Vector3D &point, const double d) {
        double X = d*point.x/-point.z;
        double Y = d*point.y/-point.z;
        return {X, Y};
}
Lines2D _3D_MatrixFunctions::doProjectionLines(const Figures3D &figures) {
        Lines2D lines;
        for (auto& f : figures) {
                for (auto& face : f.faces) {
                        img::Color color = f.color;
                        int n = face.point_indexes.size();
                        for (int i = 0; i < n; ++i) {
                                int current_index = face.point_indexes[i];
                                int next_index = face.point_indexes[(i + 1) % n];
                                Point2D pointX = doProjectionPoint(f.points[current_index], 1);
                                Point2D pointY = doProjectionPoint(f.points[next_index], 1);
                                lines.push_back(Line2D(pointX, pointY, color));
                        }
                }
        }
        return lines;
}