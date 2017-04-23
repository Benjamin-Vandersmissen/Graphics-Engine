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
                Line2D line(points[face[i]], points[face[(i+1)%size]
                ], figure.getColor());
                line.rainbow = figure.rainbow;
                if (ZBuffering){
                    line.z1 = figure[face[i]].z;
                    line.z2 = figure[face[(i+1)%size]].z;
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

unsigned int Face::operator[](unsigned int i) const{
    if (i >= this->pointIndices.size()){
        throw(std::invalid_argument("index out of range"));
    }
    else{
        return this->pointIndices[i];
    }
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

Vector3D Figure3D::operator[](unsigned int i) const{
    if (i >= this->Points.size()){
        throw(std::invalid_argument("index out of range"));
    }else{
        return this->Points[i];
    }
}

Vector3D Figure3D::getCenter() {
    Vector3D temp = Vector3D::point(0,0,0);
    for( Vector3D& point : this->Points){
        temp += point;
    }
    temp /= this->Points.size();
    return temp;
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

std::vector<Face> triangulate(const Face &face) {
    std::vector<Face> faces = {};
    for(int i = 1; i < face.getPointIndices().size()-1; i ++){
        faces.push_back(Face({face.getPointIndices()[0], face.getPointIndices()[i], face.getPointIndices()[i+1]}));
    }
    return faces;
}

void triangulate(Figure3D &figure) {
    std::vector<Face> faces = {};
    int i = 0;
    for(const Face& face : figure.getFaces()){
//        std::cerr << i << std::endl;
        std::vector<Face> newFaces = triangulate(face);
        faces.insert(faces.end(), newFaces.begin(), newFaces.end());
        i++;
    }
    figure.setFaces(faces);
}

void draw_zbuf_triangle(ZBuffer &buf, img::EasyImage &image, Vector3D &A, Vector3D &B, Vector3D &C, double d, double dx,
                        double dy, const img::Color &color) {

    Point2D newA = doProjection(A, d);
    newA.x += dx;
    newA.y += dy;

    Point2D newB = doProjection(B, d);
    newB.x += dx;
    newB.y += dy;

    Point2D newC = doProjection(C, d);
    newC.x += dx;
    newC.y += dy;

    int Ymin = round(std::min(newA.y, std::min(newB.y, newC.y))+0.5);
    int Ymax = round(std::max(newA.y, std::max(newB.y, newC.y))-0.5);

    Line2D AB = {newA,newB, img::Color()};
    Line2D AC = {newA,newC, img::Color()};
    Line2D BC = {newB,newC, img::Color()};
    Lines2D lines = {AB,AC,BC};

    img::Color color1 = color;
//    color1.blue = std::rand()%100 + 25;
//    color1.red = std::rand()%100 + 25;
//    color1.green = std::rand()%100 + 25;

    for(int y = Ymin; y <= Ymax; y++){
        double XL = std::numeric_limits<double>::infinity();
        double XR = -std::numeric_limits<double>::infinity();
        for(Line2D& l: lines){
            double X;
            try{
                X = l.getX(y);
                XL = std::min(XL, X);
                XR = std::max(XR, X);
            }catch(std::exception& e){

            }
        }
        for(int x = round(XL+0.5); x <= round(XR-0.5); x++){
            double ZGinverse = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);
            double XG = (newA.x+newB.x+newC.x)/3;
            double YG = (newA.y+newB.y+newC.y)/3;

            Vector3D U = B-A;
            Vector3D V = C-A;
            Vector3D W = Vector3D::cross(U,V);
            double k = W.x*A.x+W.y*A.y+W.z*A.z;

            double dzdx = W.x/(-d*k);
            double dzdy = W.y/(-d*k);

            double Zinverse = 1.001*ZGinverse + (x-XG)*dzdx + (y-YG)*dzdy;
            if (buf[y][x] >= Zinverse){
                buf[y][x] = Zinverse;
                image(x,y) = color1;
            }
//            else{
//                std::cerr << "HUH?" << std::endl;
//            }
        }
    }
}

bool operator==(const Face &face1, const Face &face2) {
    bool equal = face1.getPointIndices().size() == face2.getPointIndices().size();
    for(int i = 0; i < face1.getPointIndices().size(); i++){
        equal = equal && face1[i]==face2[i];
    }
    return equal;
}

