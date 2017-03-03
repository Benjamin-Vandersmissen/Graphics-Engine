//
// Created by Benjamin on 21/02/2017.
//

#ifndef ENGINE_POINT_H
#define ENGINE_POINT_H
#include <iostream>

class Point2D {
public:
    double x;
    double y;

    Point2D(double x, double y);

    Point2D();
};
std::ostream& operator<<(std::ostream& stream, Point2D& point);

#endif //ENGINE_POINT_H
