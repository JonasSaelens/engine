#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "vector3d.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>
#include <stack>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
};
using Lines2D = std::vector<Line2D>;

img::EasyImage ColorRectangle(int width, int height) {
        img::EasyImage image(width,height);
        for(unsigned int i = 0; i < width; i++)
        {
                for(unsigned int j = 0; j < height; j++)
                {
                        image(i,j).red = i;
                        image(i,j).green = j;
                        image(i,j).blue = (i+j)%256;
                }
        }
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
}
img::EasyImage Blocks(int width, int height, std::vector<double> colorWhite, std::vector<double> colorBlack, int nrXBlocks, int nrYBlocks, bool invertColors) {
        img::EasyImage image(width,height);
        int width_pixel = width/nrXBlocks;
        int height_pixel = height/nrYBlocks;
        for(unsigned int i = 0; i < width; i++)
        {
                for(unsigned int j = 0; j < height; j++)
                {
                        int Bx = i/width_pixel;
                        int By = j/height_pixel;
                        if (invertColors) {
                                if ((Bx+By)%2 == 0) {
                                        image(i,j).red = colorBlack[0]*255;
                                        image(i,j).green = colorBlack[1]*255;
                                        image(i,j).blue = colorBlack[2]*255;
                                }
                                else {
                                        image(i,j).red = colorWhite[0]*255;
                                        image(i,j).green = colorWhite[1]*255;
                                        image(i,j).blue = colorWhite[2]*255;
                                }
                        }
                        else {
                                if ((Bx+By)%2 == 0) {
                                        image(i,j).red = colorWhite[0]*255;
                                        image(i,j).green = colorWhite[1]*255;
                                        image(i,j).blue = colorWhite[2]*255;
                                }
                                else {
                                        image(i,j).red = colorBlack[0]*255;
                                        image(i,j).green = colorBlack[1]*255;
                                        image(i,j).blue = colorBlack[2]*255;
                                }
                        }
                }
        }
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
}
/*img::EasyImage Lines(int width, int height, std::string figure ,std::vector<double> backgroundcolor, std::vector<double> lineColor, int nrLines) {
        img::EasyImage image(width,height);
        if (figure =="QuarterCircle") {
                int Hs = height/(nrLines-1);
                int Ws = width/(nrLines-1);
                for(unsigned int i = 0; i < width; i++) {
                        for(unsigned int j = 0; j < height; j++)
                        {
                                if ()
                                image(i,j).red = backgroundcolor[0]*255;
                                image(i,j).green = backgroundcolor[1]*255;
                                image(i,j).blue = backgroundcolor[2]*255;
                        }
                }
        }
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
}*/

img::EasyImage drawLines2D(Lines2D &lines, const int size, std::vector<double> backgroundColor) {
        double xMin = lines[0].p1.x;
        double xMax = lines[0].p1.x;
        double yMax = lines[0].p1.y;
        double yMin = lines[0].p1.y;
        for (auto &line : lines) {
                xMin = std::min(std::min(line.p1.x,line.p2.x),xMin);
                xMax = std::max(std::max(line.p1.x,line.p2.x),xMax);
                yMax = std::max(std::max(line.p1.y,line.p2.y),yMax);
                yMin = std::min(std::min(line.p1.y,line.p2.y),yMin);
        }
        double Xrange = xMax-xMin;
        double Yrange = yMax-yMin;
        double imageX = size*(Xrange/std::max(Xrange,Yrange));
        double imageY = size*(Yrange/std::max(Xrange,Yrange));
        img::EasyImage image(lround(imageX), lround(imageY));
        for(unsigned int i = 0; i < imageX-1; i++) {
                for(unsigned int j = 0; j < imageY-1; j++)
                {
                        image(i,j).red = backgroundColor[0]*255;
                        image(i,j).green = backgroundColor[1]*255;
                        image(i,j).blue = backgroundColor[2]*255;
                }
        }
        double d = 0.95*(imageX/Xrange);
        for (auto &line : lines) {
                line.p1.x = line.p1.x * d;
                line.p1.y = line.p1.y * d;
                line.p2.x = line.p2.x * d;
                line.p2.y = line.p2.y * d;
        }
        double DCx = d*((xMin+xMax)/2);
        double DCy = d*((yMin+yMax)/2);
        double dx = (imageX/2)-DCx;
        double dy = (imageY/2)-DCy;
        for (auto &line : lines) {
                line.p1.x += dx;
                line.p1.y += dy;
                line.p2.x += dx;
                line.p2.y += dy;
                lround(line.p1.x);
                lround(line.p1.y);
                lround(line.p2.x);
                lround(line.p2.y);
                image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y, line.color);
        }
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
}

img::EasyImage LSystem2D(int size, std::vector<double> backgroundColor, std::string inputfile, std::vector<double> color){
        LParser::LSystem2D l_system;
        Lines2D list;

        std::ifstream input_stream(inputfile);
        input_stream >> l_system;
        input_stream.close();

        std::stack<Point2D> stack;
        std::stack<double> stack1;

        img::Color lineColor;
        lineColor.red = color[0]*255;
        lineColor.green = color[1]*255;
        lineColor.blue = color[2]*255;
        double radial = l_system.get_angle() * M_PI/180;
        double startingradial = l_system.get_starting_angle() * M_PI/180;
        std::set<char> alphabet = l_system.get_alphabet();
        std::vector<bool> draws;
        for (char c : alphabet) {
                bool draw = l_system.draw(c);
                draws.push_back(draw);
        }
        std::string init = l_system.get_initiator();
        std::string startingstring = init;
        std::string endstring = "";
        for (int i = 0; i < l_system.get_nr_iterations(); i++) {
                for (char k : startingstring) {
                        for (char j: alphabet) {
                                if (k == j) {
                                        endstring += l_system.get_replacement(j);
                                }
                        }
                        if (k == '+')
                                endstring += "+";
                        else if (k == '-')
                                endstring += "-";
                        else if (k == '(')
                                endstring += "(";
                        else if (k == ')')
                                endstring += ")";
                }
                startingstring = endstring;
                endstring = "";
        }
        double coorX = 0;
        double coorY = 0;
        double alpha = startingradial;
        for (char c : startingstring) {
                for (char k: alphabet) {
                        if (c==k) {
                                if (l_system.draw(k)) {
                                        double tempcoorX = coorX + cos(alpha);
                                        double tempcoorY = coorY + sin(alpha);
                                        Line2D line(Point2D(coorX, coorY), Point2D(tempcoorX,tempcoorY), lineColor);
                                        list.push_back(line);
                                        coorX = tempcoorX;
                                        coorY = tempcoorY;
                                }
                                else {
                                        coorX += cos(alpha);
                                        coorY += sin(alpha);
                                }
                        }
                }
                if (c=='+') {
                        alpha += radial;
                }
                if (c=='-') {
                        alpha -= radial;
                }
                if (c=='(') {
                        stack.push(Point2D(coorX,coorY));
                        stack1.push(alpha);
                }
                if (c==')') {
                        if (!stack.empty() && !stack1.empty()) {
                                Point2D tempPoint = stack.top();
                                alpha = stack1.top();
                                coorX = tempPoint.x;
                                coorY = tempPoint.y;
                                stack.pop();
                                stack1.pop();
                        }
                }
        }
        return drawLines2D(list, size, backgroundColor);
}

Matrix Scale(const double scale) {
        Matrix matrix;
        matrix(1,1) = scale;
        matrix(2,2) = scale;
        matrix(3,3) = scale;
        return matrix;
}
Matrix rotateX(const double angle) {
        Matrix matrix;
        double radiant = angle * (M_PI /180);
        matrix(2,2) = cos(radiant);
        matrix(2,3) = sin(radiant);
        matrix(3,2) = -sin(radiant);
        matrix(3,3) = cos(radiant);
        return matrix;
}
Matrix rotateY(const double angle) {
        Matrix matrix;
        double radiant = angle * (M_PI /180);
        matrix(1,1) = cos(radiant);
        matrix(1,3) = -sin(radiant);
        matrix(3,1) = sin(radiant);
        matrix(3,3) = cos(radiant);
        return matrix;
}
Matrix rotateZ(const double angle) {
        Matrix matrix;
        double radiant = angle * (M_PI /180);
        matrix(1,1) = cos(radiant);
        matrix(1,2) = sin(radiant);
        matrix(2,1) = -sin(radiant);
        matrix(2,2) = cos(radiant);
        return matrix;
}
Matrix translate(const Vector3D &vector) {
        Matrix matrix;
        matrix(4,1) = vector.x;
        matrix(4,2) = vector.y;
        matrix(4,3) = vector.z;
        return matrix;
}
void applyTransformation(Figure &fig, const Matrix &m) {
        for (auto &i : fig.points) {
                i *= m;
        }
}
void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
        r = sqrt(pow(point.x,2)+pow(point.y,2)+pow(point.z,2));
        theta = atan2(point.y,point.x);
        phi = acos(point.z/r);
}
Matrix eyePointTrans(const Vector3D &eyepoint) {
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
Point2D doProjectionPoint(const Vector3D &point, const double d) {
        double X = d*point.x/-point.z;
        double Y = d*point.y/-point.z;
        return Point2D(X, Y);
}
Lines2D doProjectionLines(const Figures3D &figures) {
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
int getIndex(int k, int j, int m) {
        return k * m + j;
}
img::EasyImage drawLines3D(const ini::Configuration &configuration) {
        int size = configuration["General"]["size"];
        std::vector<double> backgroundColor = configuration["General"]["backgroundcolor"];
        int nrFigures = configuration["General"]["nrFigures"];
        std::vector<double> eye = configuration["General"]["eye"];
        Figures3D figures;
        Matrix V = eyePointTrans(Vector3D::point(eye[0], eye[1], eye[2]));
        for (int i = 0; i < nrFigures; i++) {
                std::string nameFigure = "Figure" + std::to_string(i);
                std::string figureType = configuration[nameFigure]["type"];
                if (figureType == "LineDrawing") {
                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        int nrPoints = configuration[nameFigure]["nrPoints"];
                        int nrLines = configuration[nameFigure]["nrLines"];

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        Figure f;
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);
                        for (int j = 0; j < nrPoints; j++) {
                                std::string namePoint = "point" + std::to_string(j);
                                std::vector<double> point = configuration[nameFigure][namePoint];
                                f.points.push_back(Vector3D::point(point[0],point[1],point[2]));
                        }
                        for (int k = 0; k < nrLines; k++) {
                                std::string nameLine = "line" + std::to_string(k);
                                std::vector<int> line = configuration[nameFigure][nameLine];
                                Face face;
                                face.point_indexes = line;
                                f.faces.push_back(face);
                        }

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Cube") {
                        Figure f;
                        f.points.push_back(Vector3D::point(1,-1,-1));
                        f.points.push_back(Vector3D::point(-1,1,-1));
                        f.points.push_back(Vector3D::point(1,1,1));
                        f.points.push_back(Vector3D::point(-1,-1,1));
                        f.points.push_back(Vector3D::point(1,1,-1));
                        f.points.push_back(Vector3D::point(-1,-1,-1));
                        f.points.push_back(Vector3D::point(1,-1,1));
                        f.points.push_back(Vector3D::point(-1,1,1));
                        f.faces = {
                                {{0, 4, 2, 6}},
                                {{4, 1, 7, 2}},
                                {{1, 5, 3, 7}},
                                {{5, 0, 6, 3}},
                                {{6, 2, 7, 3}},
                                {{0, 5, 1, 4}}
                        };
                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Tetrahedron") {
                        Figure f;
                        f.points.push_back(Vector3D::point(1,-1,-1));
                        f.points.push_back(Vector3D::point(-1,1,-1));
                        f.points.push_back(Vector3D::point(1,1,1));
                        f.points.push_back(Vector3D::point(-1,-1,1));
                        f.points.push_back(Vector3D::point(1,1,-1));
                        f.faces = {
                                {{0, 2, 3}},
                                {{1, 3, 2}},
                                {{0, 3, 1}},
                                {{0, 2, 3}}
                        };
                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Octahedron") {
                        Figure f;
                        f.points.push_back(Vector3D::point(1,0,0));
                        f.points.push_back(Vector3D::point(0,1,0));
                        f.points.push_back(Vector3D::point(-1,0,0));
                        f.points.push_back(Vector3D::point(0,-1,0));
                        f.points.push_back(Vector3D::point(0,0,-1));
                        f.points.push_back(Vector3D::point(0,0,1));
                        f.faces = {
                                {{0, 1, 5}},
                                {{1, 2, 5}},
                                {{2, 3, 5}},
                                {{3, 0, 5}},
                                {{1, 0, 4}},
                                {{2, 1, 4}},
                                {{3, 2, 4}},
                                {{0, 3, 4}}
                        };
                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Icosahedron") {
                        Figure f;
                        f.points.push_back(Vector3D::point(0, 0, sqrt(5)/2));
                        for (int k = 2; k<7; k++) {
                                f.points.push_back(Vector3D::point(cos((k-2)*2*M_PI/5), sin((k-2)*2*M_PI/5), 0.5));
                        }
                        for (int k = 7; k<12; k++) {
                                f.points.push_back(Vector3D::point(cos(M_PI/5+(k-7)*2*M_PI/5), sin(M_PI/5+(k-7)*2*M_PI/5), -0.5));
                        }
                        f.points.push_back(Vector3D::point(0, 0, -sqrt(5)/2));
                        f.faces={
                                {{0,1,2}},{{0,2,3}},{{0,3,4}},{{0,4,5}},{{0,5,1}},
                                {{1,6,2}},{{2,6,7}},{{2,7,3}},{{3,7,8}},{{3,8,4}},
                                {{4,8,9}},{{4,9,5}},{{5,9,10}},{{5,10,1}},{{1,10,6}},
                                {{11,7,6}},{{11,8,7}},{{11,9,8}},{{11,10,9}},{{11,6,10}}};

                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Dodecahedron") {
                        Figure f1;
                        f1.points.push_back(Vector3D::point(0, 0, sqrt(5)/2));
                        for (int k = 2; k<7; k++) {
                                f1.points.push_back(Vector3D::point(cos((k-2)*2*M_PI/5), sin((k-2)*2*M_PI/5), 0.5));
                        }
                        for (int k = 7; k<12; k++) {
                                f1.points.push_back(Vector3D::point(cos(M_PI/5+(k-7)*2*M_PI/5), sin(M_PI/5+(k-7)*2*M_PI/5), -0.5));
                        }
                        f1.points.push_back(Vector3D::point(0, 0, -sqrt(5)/2));
                        f1.faces={
                                {{0,1,2}},{{0,2,3}},{{0,3,4}},{{0,4,5}},{{0,5,1}},
                                {{1,6,2}},{{2,6,7}},{{2,7,3}},{{3,7,8}},{{3,8,4}},
                                {{4,8,9}},{{4,9,5}},{{5,9,10}},{{5,10,1}},{{1,10,6}},
                                {{11,7,6}},{{11,8,7}},{{11,9,8}},{{11,10,9}},{{11,6,10}}};
                        Figure f;
                        f.faces={{{0,1,2,3,4}},{{0,5,6,7,1}},{{1,7,8,9,2}},{{2,9,10,11,3}},
                                {{3,11,12,13,4}},{{4,13,14,5,0}},{{19,18,17,16,15}},{{19,14,13,12,18}},
                                {{18,12,11,10,17}},{{17,10,9,8,16}},{{16,8,7,6,15}},{{15,6,5,14,19}}};

                        for (const Face& face : f1.faces) {
                                Vector3D p0 = f1.points[face.point_indexes[0]];
                                Vector3D p1 = f1.points[face.point_indexes[1]];
                                Vector3D p2 = f1.points[face.point_indexes[2]];

                                double resultX = (p0.x + p1.x + p2.x) / 3.0;
                                double resultY = (p0.y + p1.y + p2.y) / 3.0;
                                double resultZ = (p0.z + p1.z + p2.z) / 3.0;
                                f.points.push_back(Vector3D::point(resultX, resultY, resultZ));
                        }

                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Sphere") {
                        Figure f;
                        f.points.push_back(Vector3D::point(0, 0, sqrt(5)/2));
                        for (int k = 2; k<7; k++) {
                                f.points.push_back(Vector3D::point(cos((k-2)*2*M_PI/5), sin((k-2)*2*M_PI/5), 0.5));
                        }
                        for (int k = 7; k<12; k++) {
                                f.points.push_back(Vector3D::point(cos(M_PI/5+(k-7)*2*M_PI/5), sin(M_PI/5+(k-7)*2*M_PI/5), -0.5));
                        }
                        f.points.push_back(Vector3D::point(0, 0, -sqrt(5)/2));
                        f.faces={
                                {{0,1,2}},{{0,2,3}},{{0,3,4}},{{0,4,5}},{{0,5,1}},
                                {{1,6,2}},{{2,6,7}},{{2,7,3}},{{3,7,8}},{{3,8,4}},
                                {{4,8,9}},{{4,9,5}},{{5,9,10}},{{5,10,1}},{{1,10,6}},
                                {{11,7,6}},{{11,8,7}},{{11,9,8}},{{11,10,9}},{{11,6,10}}};

                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        int n = configuration[nameFigure]["n"];

                        for (int k=0; k < n; k++) {
                                Figure f1;
                                f1.color = f.color;
                                for (const Face& face : f.faces) {
                                        Vector3D A = f.points[face.point_indexes[0]];
                                        Vector3D B = f.points[face.point_indexes[1]];
                                        Vector3D C = f.points[face.point_indexes[2]];

                                        Vector3D D;
                                        D.x = (A.x+B.x)/2;
                                        D.y = (A.y+B.y)/2;
                                        D.z = (A.z+B.z)/2;
                                        Vector3D E;
                                        E.x = (A.x+C.x)/2;
                                        E.y = (A.y+C.y)/2;
                                        E.z = (A.z+C.z)/2;
                                        Vector3D F;
                                        F.x = (B.x+C.x)/2;
                                        F.y = (B.y+C.y)/2;
                                        F.z = (B.z+C.z)/2;

                                        int indexA = f1.points.size();
                                        f1.points.push_back(A);
                                        int indexB = f1.points.size();
                                        f1.points.push_back(B);
                                        int indexC = f1.points.size();
                                        f1.points.push_back(C);
                                        int indexD = f1.points.size();
                                        f1.points.push_back(D);
                                        int indexE = f1.points.size();
                                        f1.points.push_back(E);
                                        int indexF = f1.points.size();
                                        f1.points.push_back(F);

                                        f1.faces.push_back({{indexA, indexD, indexE}});
                                        f1.faces.push_back({{indexB, indexF, indexD}});
                                        f1.faces.push_back({{indexC, indexE, indexF}});
                                        f1.faces.push_back({{indexD, indexF, indexE}});
                                }
                                f = f1;
                        }
                        for (auto& point : f.points) {
                                point.normalise();
                        }
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Cone") {
                        Figure f;

                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        int n = configuration[nameFigure]["n"];
                        double height = configuration[nameFigure]["height"];
                        int r = 1;

                        for (int k = 0; k<n; k++) {
                                f.points.push_back(Vector3D::point(cos(2*k*M_PI/n), sin(2*k*M_PI/n), 0));
                        }
                        f.points.push_back(Vector3D::point(0,0,height));
                        std::vector<int> nConnect;
                        for (int k = 0; k<n; k ++) {
                                f.faces.push_back({{k, (k+1)%n, n}});
                                nConnect.push_back(k);
                        }
                        f.faces.push_back({nConnect});

                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Cylinder") {
                        Figure f;

                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        int n = configuration[nameFigure]["n"];
                        double height = configuration[nameFigure]["height"];
                        int r = 1;

                        for (int k = 0; k<n; k++) {
                                f.points.push_back(Vector3D::point(cos(2*k*M_PI/n), sin(2*k*M_PI/n), 0));
                        }
                        for (int k = 0; k<n; k++) {
                                f.points.push_back(Vector3D::point(cos(2*k*M_PI/n), sin(2*k*M_PI/n), height));
                        }
                        f.points.push_back(Vector3D::point(0,0,height));
                        for (int k = 0; k < n; k++) {
                                int next = (k + 1) % n;
                                f.faces.push_back({{k, next, next + n, k + n}});
                        }
                        f.faces.push_back({{n-1, 0, n,2*n-1}});

                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
                else if (figureType == "Torus") {
                        Figure f;

                        double RotateX = configuration[nameFigure]["rotateX"];
                        double RotateY = configuration[nameFigure]["rotateY"];
                        double RotateZ = configuration[nameFigure]["rotateZ"];
                        double scale = configuration[nameFigure]["scale"];
                        std::vector<double> center = configuration[nameFigure]["center"];
                        std::vector<double> color = configuration[nameFigure]["color"];
                        int n = configuration[nameFigure]["n"];
                        int m = configuration[nameFigure]["m"];
                        double r = configuration[nameFigure]["r"];
                        double R = configuration[nameFigure]["R"];

                        for (int k = 0; k < n; ++k) {
                                double theta = 2 * M_PI * k / n;
                                for (int j = 0; j < m; ++j) {
                                        double phi = 2 * M_PI * j / m;

                                        double x = (R + r * cos(phi)) * cos(theta);
                                        double y = (R + r * cos(phi)) * sin(theta);
                                        double z = r * sin(phi);

                                        f.points.push_back(Vector3D::point(x, y, z));
                                }
                        }

                        for (int k = 0; k < n; ++k) {
                                for (int j = 0; j < m; ++j) {
                                        int a = getIndex(k, j, m);
                                        int b = getIndex((k + 1) % n, j, m);
                                        int c = getIndex((k + 1) % n, (j + 1) % m, m);
                                        int d = getIndex(k, (j + 1) % m, m);
                                        f.faces.push_back({{a, b, c, d}});
                                }
                        }
                        f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                        Matrix matrix = Scale(scale) * rotateX(RotateX) * rotateY(RotateY) * rotateZ(RotateZ) * translate(Vector3D::point(center[0],center[1],center[2])) * eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                        applyTransformation(f, matrix);
                        figures.push_back(f);
                }
        }
        Lines2D list = doProjectionLines(figures);
        if (!list.empty())
                return drawLines2D(list, size, backgroundColor);
        return img::EasyImage();
}
img::EasyImage generate_image(const ini::Configuration &configuration)
{
        std::string image_type = configuration["General"]["type"];
        if (image_type == "IntroColorRectangle") {
                int image_width = configuration["ImageProperties"]["width"];
                int image_height = configuration["ImageProperties"]["height"];
                return ColorRectangle(image_width,image_height);
        }
        if (image_type == "IntroBlocks") {
                int image_width = configuration["ImageProperties"]["width"];
                int image_height = configuration["ImageProperties"]["height"];
                std::vector<double> colorWhite = configuration["BlockProperties"]["colorWhite"];
                std::vector<double> colorBlack = configuration["BlockProperties"]["colorBlack"];
                int nrXBlocks = configuration["BlockProperties"]["nrXBlocks"];
                int nrYBlocks = configuration["BlockProperties"]["nrYBlocks"];
                bool invertColors = configuration["BlockProperties"]["invertColors"];
                return Blocks(image_width,image_height,colorWhite,colorBlack,nrXBlocks,nrYBlocks,invertColors);
        }
        /*if (image_type == "IntroLines") {
                int image_width = configuration["ImageProperties"]["width"];
                int image_height = configuration["ImageProperties"]["height"];
                std::string lineType = configuration["LineProperties"]["figure"];
                std::vector<double> backgroundColor = configuration["LineProperties"]["backgroundcolor"];
                std::vector<double> lineColor = configuration["LineProperties"]["lineColor"];
                int nrLines = configuration["LineProperties"]["nrLines"];
                return Lines(image_width,image_height,lineType,backgroundColor,lineColor,nrLines);
        }*/
        if (image_type == "2DLSystem") {
                int size = configuration["General"]["size"];
                std::vector<double> backgroundColor = configuration["General"]["backgroundcolor"];
                std::string inputfile = configuration["2DLSystem"]["inputfile"];
                std::vector<double> color = configuration["2DLSystem"]["color"];
                return LSystem2D(size, backgroundColor, inputfile, color);
        }
        if (image_type == "Wireframe") {
                return drawLines3D(configuration);
        }
        return img::EasyImage();
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
/*int main()
{
        Matrix rotZ, trans;
        Vector3D p0;
        rotZ(1,1) = 0.866025404;
        rotZ(1,2) = 0.5;
        rotZ(2,1) = -0.5;
        rotZ(2,2) = rotZ(1,1);
        trans(4,1) = 1.0;
        trans(4,2) = 2.0;
        trans(4,3) = 3.0;
        Matrix combined = rotZ*trans;

        p0 = Vector3D::point(2.0, 3.0, 4.0);
        std::cout << "p0 =\t\t" << p0 << std::endl;
        Vector3D p1 = p0*rotZ;
        std::cout << "p1 = p0*rotZ =\t"
                << p1 << std::endl;
        std::cout << "p1*trans =\t"
                << (p1*trans) << std::endl;
        std::cout << "p0*combined =\t"
                << (p0*combined) << std::endl;
        return 0;
}*/

