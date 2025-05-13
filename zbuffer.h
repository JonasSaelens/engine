//
// Created by jonas on 5/11/2025.
//

#ifndef ZBUFFER_H
#define ZBUFFER_H
#include <limits>
#include <vector>


class ZBuffer : public std::vector<std::vector<double>> {
public:
    ZBuffer(const unsigned int width, const unsigned int height) {
        this->resize(width);
        for (auto& column : *this) {
            column.resize(height, std::numeric_limits<double>::infinity());
        }
    }
};



#endif //ZBUFFER_H
