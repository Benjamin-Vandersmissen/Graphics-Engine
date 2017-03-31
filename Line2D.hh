//
// Created by Benjamin on 21/02/2017.
//

#ifndef ENGINE_LINE_H
#define ENGINE_LINE_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "easy_image.hh"
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

    double z1;
    double z2;

    bool rainbow = false;

    Line2D(const Point2D &point1, const Point2D &point2, img::Color color);
    Line2D(double x1, double y1, double x2, double y2, img::Color color);

    Line2D();

    double getY(double x);
    double getX(double y);
};

typedef std::vector<Line2D> Lines2D;
img::EasyImage
draw2DLines(Lines2D &lines, const int size, const img::Color &bgColor, bool ZBuffering);
std::ostream& operator<<(std::ostream& stream, Line2D& line);

class ZBuffer: public std::vector<std::vector<double>>{
public:
    ZBuffer(const int width, const int height);
};

void draw_zbuf_line(ZBuffer &buffer, img::EasyImage &image, const unsigned int x0, const unsigned int y0,
                    const double z0, const unsigned int x1, const unsigned int y1, const double z1,
                    const img::Color &color);

void draw_zbuf_line_rainbow(ZBuffer &buffer, img::EasyImage &image, unsigned int x0, unsigned int y0, const double z0,
                            unsigned int x1, unsigned int y1, const double z1);


#endif //ENGINE_LINE_H
