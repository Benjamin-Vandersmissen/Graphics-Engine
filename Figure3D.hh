//
// Created by uauser on 3/1/17.
//

#ifndef GRAPHICS_ENGINE_FIGURE3D_HH
#define GRAPHICS_ENGINE_FIGURE3D_HH

#include <vector> //std vector
#include "vector.hh" //custom vector
#include "easy_image.hh"
#include "Face.hh"
#include "Line2D.hh"
#include <cmath>

class Figure3D {
private:
    std::vector<Face> Faces;
    std::vector<Vector3D> Points;
    img::Color color;
public:
    const std::vector<Face> &getFaces() const;

    void setFaces(const std::vector<Face> &Faces);

    const std::vector<Vector3D> &getPoints() const;

    void setPoints(const std::vector<Vector3D> &Points);

    const img::Color &getColor() const;

    void setColor(const img::Color &color);

    Figure3D(const std::vector<Face> &Faces, const std::vector<Vector3D> &Points, const img::Color &color);

    Figure3D();

    void applyTransformation(Matrix& matrix);
};
typedef std::vector<Figure3D> Figures3D;
void applyTransformation(Figures3D& figures, Matrix& matrix);
Matrix scaleFigure(const double scale);
Matrix rotateFigureX(const double angle);
Matrix rotateFigureY(const double angle);
Matrix rotateFigureZ(const double angle);
Matrix translateFigure(const Vector3D& vector);
std::ostream& operator<<(std::ostream& stream, Figure3D& figure);
Lines2D doProjection(const Figures3D& figures);
Point2D doProjection(const Vector3D& point, const double d = 1);

#endif //GRAPHICS_ENGINE_FIGURE3D_HH
