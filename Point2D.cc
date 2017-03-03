//
// Created by Benjamin on 21/02/2017.
//

#include "Point2D.hh"

Point2D::Point2D(double x, double y) : x(x), y(y) {}

Point2D::Point2D() {}

std::ostream& operator<<(std::ostream &stream, Point2D &point) {
    stream << "(" << point.x << ", " << point.y << ")";
    return stream;
}
