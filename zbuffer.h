//
// Created by jonas on 5/11/2025.
//

#ifndef ZBUFFER_H
#define ZBUFFER_H
#include <limits>
#include <vector>


class ZBuffer: public std::vector<std::vector<double> >
{
public:
    //Constructor: maakt een Z-Buffer van de correcte
    //grootte aan en initialiseert alle velden op +inf
    ZBuffer(const int width, const int height) : buffer(width, std::vector<double>(height, std::numeric_limits<double>::infinity())) {}
private:
    std::vector<std::vector<double>> buffer;
};



#endif //ZBUFFER_H
