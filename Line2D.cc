//
// Created by Benjamin on 21/02/2017.
//

#include "Line2D.hh"

//Point class
Point2D::Point2D(double x, double y) : x(x), y(y) {}

Point2D::Point2D() {}

std::ostream& operator<<(std::ostream &stream, Point2D &point) {
    stream << "(" << point.x << ", " << point.y << ")";
    return stream;
}

//Line class
Line2D::Line2D(const Point2D &point1, const Point2D &point2, img::Color color) : point1(point1), point2(point2) , color(color){}

Line2D::Line2D(double x1, double y1, double x2, double y2, img::Color color) : point1(x1, y1), point2(x2, y2), color(color){
}

Line2D::Line2D() {}

img::EasyImage
draw2DLines(Lines2D &lines, const int size, const img::Color &bgColor, bool rainbow) {
    double Xmin = lines.front().point1.x, Xmax = lines.front().point1.x, Ymin = lines.front().point1.y, Ymax = lines.front().point1.y;
    for(Line2D line: lines){
        Xmin = std::min(Xmin, std::min(line.point1.x,line.point2.x));
        Xmax = std::max(Xmax, std::max(line.point1.x,line.point2.x));
        Ymin = std::min(Ymin, std::min(line.point1.y,line.point2.y));
        Ymax = std::max(Ymax, std::max(line.point1.y,line.point2.y));
    }
    double Imagex = size * (Xmax-Xmin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
    double Imagey = size * (Ymax-Ymin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
    img::EasyImage image(roundToInt(Imagex), roundToInt(Imagey), bgColor);
    double d = 0.95 * (Imagex/ (Xmax-Xmin));
    double DCx = d * (Xmin+Xmax)/2;
    double DCy = d * (Ymin+Ymax)/2;
    double dx = (Imagex/2) - DCx;
    double dy = (Imagey/2) - DCy;
    for(Line2D line: lines){
        line.point1.x = roundToInt(d*line.point1.x+dx);
        line.point2.x = roundToInt(d*line.point2.x+dx);
        line.point1.y = roundToInt(d*line.point1.y+dy);
        line.point2.y = roundToInt(d*line.point2.y+dy);
        if (!rainbow)
        image.draw_line(line.point1.x, line.point1.y, line.point2.x, line.point2.y, line.color);
        if (rainbow)
        image.draw_line_rainbow(line.point1.x, line.point1.y, line.point2.x, line.point2.y);
    }
    return image;
}

std::ostream& operator<<(std::ostream &stream, Line2D& line) {
    stream << line.point1 << "=>" << line.point2;
    return stream;
}
