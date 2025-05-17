//
// Created by jonas on 5/17/2025.
//

#include "Triangulate.h"

std::vector<Face> Triangulate::triangulate(const Face& face) {
    std::vector<Face> triangles;
    for (int i = 1; i < face.point_indexes.size()-1; i++) {
        Face newFace = {{face.point_indexes[0],face.point_indexes[i],face.point_indexes[i+1]}};
        triangles.push_back(newFace);
    }
    return triangles;
}