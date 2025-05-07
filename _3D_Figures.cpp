//
// Created by jonas on 5/7/2025.
//

#include "_3D_Figures.h"

Figure _3D_Figures::createCube() {
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
    return f;
}
Figure _3D_Figures::createTetrahedron() {
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
    return f;
}
Figure _3D_Figures::createOctahedron() {
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
    return f;
}
Figure _3D_Figures::createIcosahedron() {
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
    return f;
}
Figure _3D_Figures::createDodecahedron() {
    Figure f1 = createIcosahedron();
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
    return f;
}
Figure _3D_Figures::createSphere(const int n) {
    Figure f = createIcosahedron();
    for (int k=0; k < n; k++) {
            Figure f1;
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
    return f;
}
Figure _3D_Figures::createCone(const int n, const double h) {
    Figure f;
    for (int k = 0; k<n; k++) {
        f.points.push_back(Vector3D::point(cos(2*k*M_PI/n), sin(2*k*M_PI/n), 0));
    }
    f.points.push_back(Vector3D::point(0,0,h));
    std::vector<int> nConnect;
    for (int k = 0; k<n; k ++) {
        f.faces.push_back({{k, (k+1)%n, n}});
        nConnect.push_back(k);
    }
    f.faces.push_back({nConnect});

    return f;
}
Figure _3D_Figures::createCylinder(const int n, const double h) {
    Figure f;
    for (int k = 0; k<n; k++) {
        f.points.push_back(Vector3D::point(cos(2*k*M_PI/n), sin(2*k*M_PI/n), 0));
    }
    for (int k = 0; k<n; k++) {
        f.points.push_back(Vector3D::point(cos(2*k*M_PI/n), sin(2*k*M_PI/n), h));
    }
    f.points.push_back(Vector3D::point(0,0,h));
    for (int k = 0; k < n; k++) {
        int next = (k + 1) % n;
        f.faces.push_back({{k, next, next + n, k + n}});
    }
    f.faces.push_back({{n-1, 0, n,2*n-1}});

    return f;
}
int getIndex(int k, int j, int m) {
    return k * m + j;
}
Figure _3D_Figures::createTorus(const double r,const double R, const int n,const int m) {
    Figure f;
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
    return f;
}
