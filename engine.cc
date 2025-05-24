#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "vector3d.h"
#include "structs.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>
#include <stack>

#include "_3D_Figures.h"
#include "_3D_MatrixFunctions.h"
#include "Triangulate.h"

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

//niet werkend
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

img::EasyImage drawLines2D(Lines2D &lines, const int size, const std::vector<double> &backgroundColor, bool zbuffering = false) {
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
        ZBuffer zbuffer(lround(imageX), lround(imageY));
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
        double dx = imageX/2-DCx;
        double dy = imageY/2-DCy;
        for (auto &line : lines) {
                line.p1.x += dx;
                line.p1.y += dy;
                line.p2.x += dx;
                line.p2.y += dy;
                if (zbuffering)
                        image.draw_zbuf_line(zbuffer,image,lround(line.p1.x), lround(line.p1.y),line.z1, lround(line.p2.x), lround(line.p2.y),line.z2, line.color);
                else
                        image.draw_line(lround(line.p1.x), lround(line.p1.y), lround(line.p2.x), lround(line.p2.y), line.color);
        }
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
}
img::EasyImage drawLines2DZbuffer(Lines2D &lines, Figures3D &f,const int size, const std::vector<double> &backgroundColor) {
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
        ZBuffer zbuffer(lround(imageX), lround(imageY));
        for(unsigned int i = 0; i < imageX-1; i++) {
                for(unsigned int j = 0; j < imageY-1; j++)
                {
                        image(i,j).red = backgroundColor[0]*255;
                        image(i,j).green = backgroundColor[1]*255;
                        image(i,j).blue = backgroundColor[2]*255;
                }
        }
        double d = 0.95*(imageX/Xrange);
        double DCx = d*((xMin+xMax)/2);
        double DCy = d*((yMin+yMax)/2);
        double dx = imageX/2-DCx;
        double dy = imageY/2-DCy;
        Triangulate t;
        for (Figure& fig : f) {
                std::vector<Face> NewFaces;
                for (Face& face : fig.faces) {
                        std::vector<Face> newFaces = t.triangulate(face);
                        NewFaces.insert(NewFaces.end(), newFaces.begin(), newFaces.end());
                }
                fig.faces = NewFaces;
                for (Face& face : fig.faces) {
                        image.draw_zbuf_triag(zbuffer,image, fig.points[face.point_indexes[0]], fig.points[face.point_indexes[1]], fig.points[face.point_indexes[2]], d, dx, dy, fig.color);
                }
        }

        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
}

img::EasyImage LSystem2D(int size, std::vector<double> backgroundColor, const std::string& inputfile, std::vector<double> color){
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
Figure LSystem3D(const std::string& inputfile){
        LParser::LSystem3D l_system;
        Figure f;

        std::ifstream input_stream(inputfile);
        input_stream >> l_system;
        input_stream.close();

        std::stack<Vector3D> stack;
        Vector3D H = Vector3D::vector(1,0,0);
        Vector3D L = Vector3D::vector(0,1,0);
        Vector3D U = Vector3D::vector(0,0,1);
        std::stack<Vector3D> stack1;
        std::stack<Vector3D> stack2;
        std::stack<Vector3D> stack3;

        double radial = l_system.get_angle() * M_PI/180;
        std::set<char> alphabet = l_system.get_alphabet();
        std::string init = l_system.get_initiator();
        std::string startingstring = init;
        std::string endstring;
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
                        else if (k == '^')
                                endstring += "^";
                        else if (k == '&')
                                endstring += "&";
                        else if (k == '\\')
                                endstring += "\\";
                        else if (k == '/')
                                endstring += "/";
                        else if (k == '|')
                                endstring += "|";
                }
                startingstring = endstring;
                endstring = "";
        }
        //std::cout << startingstring << std::endl;
        Vector3D coor = Vector3D::point(0,0,0);
        f.points.push_back(coor);
        int teller = 1;
        for (char c : startingstring) {
                Vector3D copyH = H;
                Vector3D copyL = L;
                Vector3D copyU = U;
                for (char k: alphabet) {
                        if (c==k) {
                                if (l_system.draw(k)) {
                                        coor += H;
                                        f.points.push_back(coor);
                                        f.faces.push_back({{teller, teller-1}});
                                        teller++;
                                }
                                else {
                                        coor += H;
                                        f.points.push_back(coor);
                                        teller++;
                                }
                        }
                }
                if (c == '+') {
                        H = H*cos(radial) + L*sin(radial);
                        L = -copyH*sin(radial) + L*cos(radial);
                }
                else if (c == '-') {
                        H = H*cos(-radial) + L*sin(-radial);
                        L = -copyH*sin(-radial) + L*cos(-radial);
                }
                else if (c == '^') {
                        H = H*cos(radial) + U*sin(radial);
                        U = -copyH*sin(radial) + U*cos(radial);
                }
                else if (c == '&') {
                        H = H*cos(-radial) + U*sin(-radial);
                        U = -copyH*sin(-radial) + U*cos(-radial);
                }
                else if (c == '\\') {
                        L = L*cos(radial) - U*sin(radial);
                        U = copyL*sin(radial) + U*cos(radial);
                }
                else if (c == '/') {
                        L = L*cos(-radial) - U*sin(-radial);
                        U = copyL*sin(-radial) + U*cos(-radial);
                }
                else if (c == '|') {
                        H = -H;
                        L = -L;
                }
                else if (c == '(') {
                        stack.push(coor);
                        stack1.push(H);
                        stack2.push(L);
                        stack3.push(U);
                }
                else if (c == ')') {
                        if (!stack.empty() && !stack1.empty() && !stack2.empty() && !stack3.empty()) {
                                coor = stack.top();
                                H = stack1.top();
                                L = stack2.top();
                                U = stack3.top();
                                stack.pop();
                                stack1.pop();
                                stack2.pop();
                                stack3.pop();
                                f.points.push_back(coor);
                                teller++;
                        }
                }
        }
        return f;
}

Figures3D generateFractal(Figure& fig,const int nr_iterations, const double scale) {
        std::vector<Figure> figures = {fig};
        for (int i = 0; i<nr_iterations; i++) {
                std::vector<Figure> currentNewFigures;
                for (Figure& current : figures) {
                        int index = 0;
                        for (auto& point : current.points) {
                                Figure temp = current;
                                _3D_MatrixFunctions::applyTransformation(temp, _3D_MatrixFunctions::Scale(1/scale));
                                _3D_MatrixFunctions::applyTransformation(temp, _3D_MatrixFunctions::translate(point - temp.points[index]));
                                index++;
                                currentNewFigures.push_back(temp);
                        }
                }
                figures = currentNewFigures;
        }

        Figures3D fractal;
        for (Figure& figure : figures) {
                fractal.push_back(figure);
        }
        return fractal;
}

img::EasyImage drawLines3D(const ini::Configuration &configuration, bool zbuffering = false, bool triangle = false) {
        int size = configuration["General"]["size"];
        std::vector<double> backgroundColor = configuration["General"]["backgroundcolor"];
        int nrFigures = configuration["General"]["nrFigures"];
        std::vector<double> eye = configuration["General"]["eye"];
        Figures3D figures;
        bool containsFractal = false;
        std::vector<bool> fractal(nrFigures, false);
        std::vector<int> nrIt(nrFigures, 0);
        std::vector<double> fractalScale(nrFigures, 0);
        for (int i = 0; i < nrFigures; i++) {
                Figure f;
                std::string nameFigure = "Figure" + std::to_string(i);
                std::string figureType = configuration[nameFigure]["type"];
                double RotateX = configuration[nameFigure]["rotateX"];
                double RotateY = configuration[nameFigure]["rotateY"];
                double RotateZ = configuration[nameFigure]["rotateZ"];
                double scale = configuration[nameFigure]["scale"];
                std::vector<double> center = configuration[nameFigure]["center"];
                std::vector<double> color = configuration[nameFigure]["color"];
                if (figureType == "FractalCube" || figureType == "FractalTetrahedron" || figureType == "FractalOctahedron" || figureType == "FractalIcosahedron" || figureType == "FractalDodecahedron" || figureType == "FractalLineDrawing" || figureType == "FractalSphere" || figureType == "FractalCone" || figureType == "FractalCylinder" || figureType == "FractalTorus" || figureType == "Fractal3DLSystem") {
                        containsFractal = true;
                        fractal[i] = true;
                        int nrItTemp = configuration[nameFigure]["nrIterations"];
                        double fractalScaleTemp = configuration[nameFigure]["fractalScale"];
                        nrIt[i] = nrItTemp;
                        fractalScale[i] = fractalScaleTemp;
                }
                if (figureType == "LineDrawing" || figureType == "FractalLineDrawing") {
                        int nrPoints = configuration[nameFigure]["nrPoints"];
                        int nrLines = configuration[nameFigure]["nrLines"];
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
                }
                else if (figureType == "Cube" || figureType == "FractalCube") {
                        f = _3D_Figures::createCube();
                }
                else if (figureType == "Tetrahedron" || figureType == "FractalTetrahedron") {
                        f = _3D_Figures::createTetrahedron();
                }
                else if (figureType == "Octahedron" || figureType == "FractalOctahedron") {
                        f = _3D_Figures::createOctahedron();
                }
                else if (figureType == "Icosahedron" || figureType == "FractalIcosahedron") {
                        f = _3D_Figures::createIcosahedron();
                }
                else if (figureType == "Dodecahedron" || figureType == "FractalDodecahedron") {
                        f = _3D_Figures::createDodecahedron();
                }
                else if (figureType == "Sphere" || figureType == "FractalSphere") {
                        int n = configuration[nameFigure]["n"];

                        f = _3D_Figures::createSphere(n);
                }
                else if (figureType == "Cone" || figureType == "FractalCone") {
                        int n = configuration[nameFigure]["n"];
                        double height = configuration[nameFigure]["height"];

                        f = _3D_Figures::createCone(n, height);

                }
                else if (figureType == "Cylinder" || figureType == "FractalCylinder") {
                        int n = configuration[nameFigure]["n"];
                        double height = configuration[nameFigure]["height"];

                        f = _3D_Figures::createCylinder(n, height);

                }
                else if (figureType == "Torus" || figureType == "FractalTorus") {
                        int n = configuration[nameFigure]["n"];
                        int m = configuration[nameFigure]["m"];
                        double r = configuration[nameFigure]["r"];
                        double R = configuration[nameFigure]["R"];

                        f = _3D_Figures::createTorus(r,R,n,m);
                }
                else if (figureType == "3DLSystem" || figureType == "Fractal3DLSystem") {
                        std::string inputFile = configuration[nameFigure]["inputfile"];
                        f = LSystem3D(inputFile);
                }
                f.color = img::Color(color[0]*255, color[1]*255, color[2]*255);

                Matrix matrix = _3D_MatrixFunctions::Scale(scale) * _3D_MatrixFunctions::rotateX(RotateX) * _3D_MatrixFunctions::rotateY(RotateY) * _3D_MatrixFunctions::rotateZ(RotateZ) * _3D_MatrixFunctions::translate(Vector3D::point(center[0],center[1],center[2])) * _3D_MatrixFunctions::eyePointTrans(Vector3D::point(eye[0],eye[1],eye[2]));;

                _3D_MatrixFunctions::applyTransformation(f, matrix);
                figures.push_back(f);
        }
        Figures3D finalFigures;

        if (containsFractal) {
                int index = 0;
                for (auto& fig : figures) {
                        if (fractal[index]) {
                                Figures3D temp = generateFractal(fig, nrIt[index], fractalScale[index]);
                                for (Figure& temp1 : temp) {
                                        finalFigures.push_back(temp1);
                                }
                        }
                        else {
                                finalFigures.push_back(fig);
                        }
                        index++;
                }
        }
        if (!finalFigures.empty()) {
                figures = finalFigures;
        }
        Lines2D list = _3D_MatrixFunctions::doProjectionLines(figures);
        if (!list.empty()) {
                if (zbuffering) {
                        if (triangle)
                                return drawLines2DZbuffer(list, figures, size, backgroundColor);
                        return drawLines2D(list, size, backgroundColor, true);
                }
                return drawLines2D(list, size, backgroundColor);
        }
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
        if (image_type == "ZBufferedWireframe") {
                return drawLines3D(configuration, true);
        }
        if (image_type == "ZBuffering") {
                return drawLines3D(configuration, true, true);
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