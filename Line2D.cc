//
// Created by Benjamin on 21/02/2017.
//

#include <limits>
#include <cassert>
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
draw2DLines(Lines2D &lines, const int size, const img::Color &bgColor, bool rainbow, bool ZBuffering) {
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
    ZBuffer buffer(roundToInt(Imagex), roundToInt(Imagey));

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
        if (ZBuffering){
            draw_zbuf_line(buffer, image, line.point1.x, line.point1.y, line.z1, line.point2.x , line.point2.y, line.z2, line.color);
        }
        else if (!rainbow )
        image.draw_line(line.point1.x, line.point1.y, line.point2.x, line.point2.y, line.color);
        else
        image.draw_line_rainbow(line.point1.x, line.point1.y, line.point2.x, line.point2.y);
    }
    return image;
}

std::ostream& operator<<(std::ostream &stream, Line2D& line) {
    stream << line.point1 << "=>" << line.point2;
    return stream;
}

ZBuffer::ZBuffer(const int width, const int height) {
    const double infinity = std::numeric_limits<double>::infinity();
    for(int i = 0; i < height; i++){
        this->push_back({});
        for (int j = 0; j < width; j++){
            this->at(i).push_back(infinity);
        }
    }
}

void draw_zbuf_line(ZBuffer &buffer, img::EasyImage &image,  unsigned int x0,  unsigned int y0,
                    const double z0,  unsigned int x1,  unsigned int y1, const double z1,
                    const img::Color &color) {
    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());


    if (x0 == x1)
    {
        unsigned int ymin = std::min(y0,y1);
        unsigned int ymax = std::max(y0,y1);
        double zmin = (ymin == y0) ? z0 : z1;
        double zmax = (ymax == y0) ? z0 : z1;
        //special case for x0 == x1
        for (unsigned int i = ymin; i <= ymax; i++)
        {
            double z = zmin + (i-ymin)*(zmax-zmin)/(ymax-ymin);
            if (buffer[i][x0] > 1/z) {
                buffer[i][x0] = 1/z;
                image(x0, i) = color;
            }
            else {
                //std::cerr << buffer[i][x0] << ' ' << 1/z << std::endl;
            }
        }
    }
    else if (y0 == y1)
    {
        unsigned int xmin = std::min(x0,x1);
        unsigned int xmax = std::max(x0,x1);
        double zmin = (xmin == x0) ? z0 : z1;
        double zmax = (xmax == x0) ? z0 : z1;
        //special case for y0 == y1
        for (unsigned int i = xmin; i <= xmax; i++)
        {
            double z = zmin + (i-xmin)*(zmax-zmin)/(xmax-xmin);
            if (buffer[y0][i] > 1/z) {
                buffer[y0][i] = 1/z;
                image(i, y0) = color;
            }
            else{
            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double z = z0 + i*(z1-z0)/(x1-x0);
                if (buffer[(unsigned int) round(y0 + m * i)][x0+i] > 1/z ) {
                    buffer[(unsigned int) round(y0 + m * i)][x0+i] = 1/z;
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }else{
                    //std::cerr << buffer[x0+i][(unsigned int) round(y0 + m * i)] << ' ' << 1/z << std::endl;
                }
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double z = z0 + i*(z1-z0)/(y1-y0);
                if(buffer[y0+i][(unsigned int) round(x0 +(i / m))] > 1/z) {
                    buffer[y0+i][(unsigned int) round(x0 +(i / m))] = 1/z;
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }else{
//                    std::cerr << buffer[y0+i][(unsigned int) round(x0 +(i / m))] << ' ' << 1/z << std::endl;

                }
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double z = z0 + i*(z1-z0)/(y0-y1);
                if (buffer[y0-i][(unsigned int) round(x0 - (i / m))] > 1/z){
                    buffer[y0-i][(unsigned int) round(x0 - (i / m))] = 1/z;
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }else{
//                    std::cerr << buffer[y0-i][(unsigned int) round(x0 - (i / m))] << ' ' << 1/z << std::endl;
                }
            }
        }
    }
}
