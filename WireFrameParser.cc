//
// Created by Benjamin on 01/03/2017.
//

#include "WireFrameParser.hh"

WireFrameParser::WireFrameParser(const ini::Configuration &configuration) {
    int size = configuration["General"]["size"].as_int_or_die();
    img::Color bgcolor = img::Color(extractColor(configuration["General"]["backgroundcolor"].as_double_tuple_or_die()));
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    std::vector<double> EyeCoords = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eye;
    eye = Vector3D::point(EyeCoords[0], EyeCoords[1], EyeCoords[2]);
    Figures3D figures;
    for(int i = 0; i < nrFigures; i++){
        std::string name = "Figure" + std::to_string(i);
        std::string type = configuration[name.c_str()]["type"].as_string_or_die();
        double rotateX = toRadial(configuration[name.c_str()]["rotateX"].as_double_or_die());
        double rotateY = toRadial(configuration[name.c_str()]["rotateY"].as_double_or_die());
        double rotateZ = toRadial(configuration[name.c_str()]["rotateZ"].as_double_or_die());
        double scale = configuration[name.c_str()]["scale"].as_double_or_die();
        std::vector<double> CenterCoords = configuration[name.c_str()]["center"].as_double_tuple_or_die();
        Vector3D CenterPoint;
        CenterPoint = Vector3D::point(CenterCoords[0], CenterCoords[1], CenterCoords[2]);
        img::Color color = img::Color(extractColor(configuration[name.c_str()]["color"].as_double_tuple_or_die()));

        Matrix m;

        Vector3D vec;
        if (type == "LineDrawing"){
            figures.push_back(this->parseLinedrawing(configuration, name, color));
        }
        else if(type == "Cube"){
            figures.push_back(this->parseCube(color));
        }
        else if(type == "Tetrahedron"){
            figures.push_back(this->parseTetrahedron(color));
        }
        else if (type == "Octahedron"){
            figures.push_back(this->parseOctahedron(color));
        }
        else if (type == "Icosahedron"){
            figures.push_back(this->parseIcosahedron(color));
        }

        m = scaleFigure(scale);
        figures.back().applyTransformation(m); //apply scaling

        m = rotateFigureX(rotateX);
        figures.back().applyTransformation(m); //rotate X-axis

        m = rotateFigureY(rotateY);
        figures.back().applyTransformation(m); //rotate Y-axis

        m = rotateFigureZ(rotateZ);
        figures.back().applyTransformation(m); //rotate Z-axis

        vec = Vector3D::vector(+CenterPoint.x, +CenterPoint.y, +CenterPoint.z); //translate the whole Figure if the centerpoint isn't (0,0,0);
        m = translateFigure(vec);
        figures.back().applyTransformation(m);
        CenterPoint = CenterPoint.point(0,0,0); //figure is translated so centerpoint is set to (0,0,0)
    }
    double r = std::sqrt((pow(eye.x,2) + pow(eye.y,2) + pow(eye.z,2)));
    double theta = std::atan2(eye.y,eye.x);
    double phi = std::acos(eye.z/r);
    toPolar(eye, r, theta, phi);

    Matrix eyeMatrix;
    eyeMatrix(1,1) = -sin(theta);
    eyeMatrix(1,2) = -cos(theta) * cos(phi);
    eyeMatrix(1,3) = cos(theta) * sin(phi);

    eyeMatrix(2,1) = cos(theta);
    eyeMatrix(2,2) = -sin(theta) * cos(phi);
    eyeMatrix(2,3) = sin(theta)*sin(phi);

    eyeMatrix(3,2) = sin(phi);
    eyeMatrix(3,3) = cos(phi);

    eyeMatrix(4,3) = -r;

    applyTransformation(figures, eyeMatrix);
    Lines2D lines = doProjection(figures);
    img::EasyImage image;
    image = draw2DLines(lines, size, bgcolor);
    this->image = image;
}

Figure3D WireFrameParser::parseLinedrawing(const ini::Configuration &configuration, std::string &name,
                                           img::Color &color) {
    int nrPoints = configuration[name.c_str()]["nrPoints"].as_int_or_die();
    int nrLines = configuration[name.c_str()]["nrLines"].as_int_or_die();
    std::string str;
    std::vector<Vector3D> Points;
    for(int i = 0; i < nrPoints; i++){
        str = "point" + std::to_string(i);
        std::vector<double> Coordinates = configuration[name.c_str()][str.c_str()].as_double_tuple_or_die();
        Vector3D point;
        point = Vector3D::point(Coordinates[0], Coordinates[1], Coordinates[2]);
        Points.push_back(point);
    }
    std::vector<Face> Faces;
    for(int i = 0; i < nrLines; i++){
        str = "line" + std::to_string(i);
        std::vector<int> indices = configuration[name.c_str()][str.c_str()].as_int_tuple_or_die();
        for(int j = 0; j < indices.size(); j++){
            if (indices[j] >= Points.size() || indices[j] < 0){
                std::cerr << "Index out of bounds" << std::endl;
                throw(std::invalid_argument("Index out of bounds"));
            }
        }
        Faces.push_back(Face(indices));
    }
    Figure3D figure(Faces, Points, color);
    return figure;
}

const img::EasyImage &WireFrameParser::getImage() const {
    return image;
}

Figure3D WireFrameParser::parseCube(img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points= {Vector3D::point(1,-1,-1), Vector3D::point(-1,1,-1), Vector3D::point(1,1,1), Vector3D::point(-1,-1,1),
                                   Vector3D::point(1,1,-1), Vector3D::point(-1,-1,-1), Vector3D::point(1,-1,1), Vector3D::point(-1,1,1)};
    figure.setPoints(points);
    figure.setColor(color);
    std::vector<Face> faces = {Face({1,5,3,7}),Face({5,2,8,3}), Face({2,6,4,8}), Face({6,1,7,4}), Face({7,3,8,4}), Face({1,6,2,5})};
    figure.setFaces(faces);
    return figure;
}

Figure3D WireFrameParser::parseTetrahedron(img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(1,-1,-1), Vector3D::point(-1,1,-1), Vector3D::point(1,1,1), Vector3D::point(-1,-1,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({1,2,3}), Face({2,4,3}), Face({1,4,2}), Face({1,3,4})};
    figure.setFaces(faces);
    figure.setColor(color);
    return figure;
}

Figure3D WireFrameParser::parseOctahedron(img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(1,0,0), Vector3D::point(0,1,0), Vector3D::point(-1,0,0), Vector3D::point(0,-1,0),
                                    Vector3D::point(0,0,-1), Vector3D::point(0,0,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({1,2,6}), Face({2,3,6}), Face({3,4,6}), Face({4,1,6}),
                               Face({2,1,5}), Face({3,2,5}), Face({4,3,5}), Face({1,4,5})};
    figure.setFaces(faces);
    figure.setColor(color);
    return figure;
}

Figure3D WireFrameParser::parseIcosahedron(img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(0,0, std::sqrt(5)/2)};
    for(int i = 2; i < 12; i ++){
        Vector3D point;
        if(i <= 6){
            point = point.point(cos((i-2)*2*M_PI/5), sin((i-2)*2*M_PI/5), 0.5);
        }
        else{
            point = point.point(cos(M_PI/5 + (i-7)*2*M_PI/5), sin(M_PI/5 +(i-7)*2*M_PI/5), -0.5);
        }
        points.push_back(point);
    }
    points.push_back(Vector3D::point(0,0,-sqrt(5)/2));
    figure.setPoints(points);
    std::vector<Face> faces = {Face({1,2,3}), Face({1,3,4}), Face({1,4,5}), Face({1,5,6}), Face({1,6,2}),
                               Face({2,7,3}), Face({3,7,8}), Face({3,8,4}), Face({4,8,9}), Face({4,9,5}),
                               Face({5,9,10}), Face({5,10,6}), Face({6,10,11}), Face({6,11,2}), Face({2,11,7}),
                               Face({12,8,7}), Face({12,9,8}), Face({12,10,9}), Face({12,11,10}), Face({12,7,11})};
    figure.setFaces(faces);
    figure.setColor(color);
    return figure;
}
