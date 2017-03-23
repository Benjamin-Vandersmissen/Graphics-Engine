//
// Created by uauser on 3/1/17.
//

#include "Figure3D.hh"

//transformation functions
Matrix scaleFigure(const double scale) {
    Matrix matrix;
    matrix(1,1) = scale;
    matrix(2,2) = scale;
    matrix(3,3) = scale;
    return matrix;
}

Matrix rotateFigureX(const double angle) {
    Matrix matrix;
    matrix(2,2) = std::cos(angle);
    matrix(2,3) = std::sin(angle);
    matrix(3,2) = -std::sin(angle);
    matrix(3,3) = std::cos(angle);
    return matrix;
}

Matrix rotateFigureY(const double angle) {
    Matrix matrix;
    matrix(1,1) = std::cos(angle);
    matrix(1,3) = -std::sin(angle);
    matrix(3,1) = std::sin(angle);
    matrix(3,3) = std::cos(angle);
    return matrix;
}

Matrix rotateFigureZ(const double angle) {
    Matrix matrix;
    matrix(1,1) = std::cos(angle);
    matrix(1,2) = std::sin(angle);
    matrix(2,1) = -std::sin(angle);
    matrix(2,2) = std::cos(angle);
    return matrix;
}

Matrix translateFigure(const Vector3D &vector) {
    Matrix matrix;
    matrix(4,1) = vector.x;
    matrix(4,2) = vector.y;
    matrix(4,3) = vector.z;
    return matrix;
}

void applyTransformation(Figures3D &figures, Matrix &matrix) {
    for(Figure3D& figure : figures){
        figure.applyTransformation(matrix);
    }
}

void Figure3D::applyTransformation(Matrix &matrix) {
    for(Vector3D& point : this->Points){
        point *= matrix;
    }
}


//projection functions
Lines2D doProjection(const Figures3D &figures, bool ZBuffering) {
    Lines2D lines;
    for(const Figure3D& figure : figures){
        std::vector<Point2D> points;
        for(const Vector3D& point : figure.getPoints()){
            Point2D newpoint = doProjection(point);
            points.push_back(newpoint);
        }
        for(const Face face: figure.getFaces()){
            for(int i = 0; i < face.getPointIndices().size(); i++){
                //connect all points from the face
                int size = face.getPointIndices().size();
                Line2D line(points[face.getPointIndices()[i]], points[face.getPointIndices()[(i+1)%size]
                ], figure.getColor());
                line.rainbow = figure.rainbow;
                if (ZBuffering){
                    line.z1 = figure.getPoints()[face.getPointIndices()[i]-1].z;
                    line.z2 = figure.getPoints()[face.getPointIndices()[(i+1)%size]-1].z;
                }
                lines.push_back(line);
                if (size == 2){ //size = 2 -> normally : draw_line(p0, p1) and draw_line(p1,p0) => avoid this behaviour
                    break;
                }
            }
        }
    }
    return lines;
}

Point2D doProjection(const Vector3D &point, const double d) {
    return Point2D(d*point.x/-point.z, d*point.y/-point.z);
}

// Face class
Face::Face(const std::vector<int> &pointIndices) : pointIndices(pointIndices) {}

Face::Face() {}

const std::vector<int> &Face::getPointIndices() const {
    return pointIndices;
}

//Figure3D class
Figure3D::Figure3D(const std::vector<Face> &Faces, const std::vector<Vector3D> &Points, const img::Color &color)
        : Faces(Faces), Points(Points), color(color) {}


const std::vector<Face> &Figure3D::getFaces() const {
    return Faces;
}

void Figure3D::setFaces(const std::vector<Face> &Faces) {
    Figure3D::Faces = Faces;
}

const std::vector<Vector3D> &Figure3D::getPoints() const {
    return Points;
}

void Figure3D::setPoints(const std::vector<Vector3D> &Points) {
    Figure3D::Points = Points;
}

const img::Color &Figure3D::getColor() const {
    return color;
}

void Figure3D::setColor(const img::Color &color) {
    Figure3D::color = color;
}

Figure3D::Figure3D() {}

Vector3D Figure3D::getCenter(int face) {
    if (face < 0 || face > this->Faces.size()){
        throw(std::invalid_argument("index out of bounds"));
    }
    Vector3D point = Vector3D::point(0,0,0);
    for(int index: this->Faces[face].getPointIndices()){
        point += this->Points[index];
    }
    point /= this->Faces[face].getPointIndices().size();
    return point;
}


std::ostream &operator<<(std::ostream &stream, const Figure3D& figure) {
    for(Face face : figure.getFaces()) {
        for (int i = 0; i < face.getPointIndices().size(); i++) {
            stream << figure.getPoints()[face.getPointIndices()[i]-1];
            stream << "->";
            stream << figure.getPoints()[face.getPointIndices()[(i+1)%face.getPointIndices().size()]-1] << std::endl;
        }
    }
    return stream;
}
