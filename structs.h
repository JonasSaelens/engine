//
// Created by jonas on 5/7/2025.
//

#ifndef STRUCTS_H
#define STRUCTS_H

#include "easy_image.h"
#include "vector3d.h"

#include <list>

class Face
{
public:
    //De indexen refereren naar
    //punten in de ‘points’ vector
    //van de Figure-klasse
    std::vector<int> point_indexes;
};
class Figure
{
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
};
typedef std::list<Figure> Figures3D;

class Point2D {
public:
    Point2D(double x, double y)
            : x(x),
              y(y) {
    }

    double x;
    double y;
};
class Line2D {
public:
    Line2D(const Point2D &p1, const Point2D &p2, const img::Color &color)
            : p1(p1),
              p2(p2),
              color(color) {
    }

    Point2D p1;
    Point2D p2;
    img::Color color;

    double z1;
    double z2;
};
using Lines2D = std::vector<Line2D>;



#endif //STRUCTS_H
