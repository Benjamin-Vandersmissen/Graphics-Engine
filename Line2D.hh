//
// Created by Benjamin on 21/02/2017.
//

#ifndef ENGINE_LINE_H
#define ENGINE_LINE_H

#include <vector>
#include "easy_image.hh"
#include <cmath>
#include "UsefulFunctions.hh"

class Point2D {
public:
    double x;
    double y;

    Point2D(double x, double y);

    Point2D();
};
std::ostream& operator<<(std::ostream& stream, Point2D& point);

class Line2D {
public:
    Point2D point1;
    Point2D point2;
    img::Color color;

    Line2D(const Point2D &point1, const Point2D &point2, img::Color color);
    Line2D(double x1, double y1, double x2, double y2, img::Color color);

    Line2D();
};

typedef std::vector<Line2D> Lines2D;
img::EasyImage
draw2DLines(Lines2D &lines, const int size, const img::Color &bgColor, bool rainbow = false);
std::ostream& operator<<(std::ostream& stream, Line2D& line);


#endif //ENGINE_LINE_H
