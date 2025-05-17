//
// Created by jonas on 5/17/2025.
//

#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include "structs.h"

class Triangulate {
public:
    std::vector<Face> triangulate(const Face& face);
};



#endif //TRIANGULATE_H
