//
// Created by Benjamin on 01/03/2017.
//

#include "WireFrameParser.hh"

WireFrameParser::WireFrameParser(const ini::Configuration &configuration, unsigned int ZBuffering) : ZBuffering(ZBuffering) {
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

        bool rainbow = configuration[name.c_str()]["rainbow"].as_bool_or_default(false);

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
        else if (type == "Dodecahedron"){
            figures.push_back(this->parseDodecahedron(color));
        }
        else if (type == "Cone"){
            figures.push_back(this->parseCone(configuration, name, color));
        }
        else if (type == "Cylinder"){
            figures.push_back(this->parseCylinder(configuration, name, color));
        }
        else if (type == "Sphere"){
            figures.push_back(this->parseSphere(configuration, name, color));
        }
        else if (type == "Torus"){
            figures.push_back(this->parseTorus(configuration, name, color));
        }
//        else if (type == "Mobius"){
//            figures.push_back(this->parseMobius(configuration, name, color));
//        }
        else if (type == "3DLSystem"){
            figures.push_back(this->parse3DLsystem(configuration, name, color));
        }
//        else if (type == "NavelTorus"){
//            figures.push_back(this->parseNavelTorus(configuration, name, color));
//        }
        figures.back().rainbow = rainbow;

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
    Lines2D lines = doProjection(figures, ZBuffering != 0);
    img::EasyImage image;
    if (ZBuffering < 2)
        image = draw2DLines(lines, size, bgcolor, ZBuffering==1);
    else if (ZBuffering == 2){
        double Xmin = lines.front().point1.x, Xmax = lines.front().point1.x, Ymin = lines.front().point1.y, Ymax = lines.front().point1.y;
        for(Line2D line: lines){
            Xmin = std::min(Xmin, std::min(line.point1.x,line.point2.x));
            Xmax = std::max(Xmax, std::max(line.point1.x,line.point2.x));
            Ymin = std::min(Ymin, std::min(line.point1.y,line.point2.y));
            Ymax = std::max(Ymax, std::max(line.point1.y,line.point2.y));
        }
        double Imagex = size * (Xmax-Xmin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
        double Imagey = size * (Ymax-Ymin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
        image = img::EasyImage(roundToInt(Imagex), roundToInt(Imagey), bgcolor);
        ZBuffer buffer(roundToInt(Imagex), roundToInt(Imagey));
        double d = 0.95 * (Imagex/ (Xmax-Xmin));
        double DCx = d * (Xmin+Xmax)/2;
        double DCy = d * (Ymin+Ymax)/2;
        double dx = (Imagex/2) - DCx;
        double dy = (Imagey/2) - DCy;
        for (Figure3D& figure: figures){
            triangulate(figure);
//            std::cerr << figure.getFaces().size();
            for(const Face& face : figure.getFaces()){
                Vector3D pointA = figure[face[0]];
                Vector3D pointB = figure[face[1]];
                Vector3D pointC = figure[face[2]];
                draw_zbuf_triangle(buffer, image, pointA, pointB, pointC, d, dx, dy, figure.getColor());
            }
        }
    }
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
    std::vector<Face> faces = {Face({0,4,2,6}),Face({4,1,7,2}), Face({1,5,3,7}), Face({5,0,6,3}), Face({6,2,7,3}), Face({0,5,1,4})};
    figure.setFaces(faces);
    return figure;
}

Figure3D WireFrameParser::parseTetrahedron(img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(1,-1,-1), Vector3D::point(-1,1,-1), Vector3D::point(1,1,1), Vector3D::point(-1,-1,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({0,1,2}), Face({1,3,2}), Face({0,3,1}), Face({0,2,3})};
    figure.setFaces(faces);
    figure.setColor(color);
    return figure;
}

Figure3D WireFrameParser::parseOctahedron(img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(1,0,0), Vector3D::point(0,1,0), Vector3D::point(-1,0,0), Vector3D::point(0,-1,0),
                                    Vector3D::point(0,0,-1), Vector3D::point(0,0,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({0,1,5}), Face({1,2,5}), Face({2,3,5}), Face({3,0,5}),
                               Face({2,1,5}), Face({2,1,4}), Face({3,2,4}), Face({0,3,4})};
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
    std::vector<Face> faces = {Face({0,1,2}), Face({0,2,3}), Face({0,3,4}), Face({0,4,5}), Face({0,5,1}),
                               Face({1,6,2}), Face({2,6,7}), Face({2,7,3}), Face({3,7,8}), Face({3,8,4}),
                               Face({4,8,9}), Face({4,9,5}), Face({5,9,10}), Face({5,10,1}), Face({1,10,6}),
                               Face({11,7,6}), Face({11,8,7}), Face({11,9,8}), Face({11,10,9}), Face({11,6,10})};
    figure.setFaces(faces);
    figure.setColor(color);
    return figure;
}

Figure3D WireFrameParser::parseDodecahedron(img::Color &color) {
    Figure3D figure = parseIcosahedron(color);
    Figure3D figure2;
    std::vector<Vector3D> points = {};
    for(int i = 0; i < figure.getFaces().size();i++){
        points.push_back(figure.getCenter(i));
    }
    figure2.setPoints(points);
    std::vector<Face> faces = {Face({0,1,2,3,4}), Face({0,5,6,7,1}), Face({1,7,8,9,2}), Face({2,9,10,11,3}),
                               Face({3,11,12,13,4}), Face({4,13,14,5,0}), Face({19,18,17,16,15}),
                               Face({19,14,13,12,18}), Face({18,12,11,10,17}), Face({17,10,9,8,16}),
                               Face({16,8,7,6,15}), Face({15,6,5,14,19})};

    figure2.setFaces(faces);
    figure2.setColor(color);
    return figure2;
}

Figure3D WireFrameParser::parseCone(const ini::Configuration &configuration, std::string &name, img::Color &color) {
    double height = configuration[name]["height"].as_double_or_die();
    int n = configuration[name]["n"].as_int_or_die();
    Figure3D figure;
    Vector3D top = Vector3D::point(0,0,height);
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    std::vector<int> circleIndices = {};
    for(int i = 0; i < n; i++){
        points.push_back(Vector3D::point(cos(2*i*M_PI/n), sin(2*i*M_PI/n) , 0));
        faces.push_back(Face({i, (i+1)%n, n}));
        circleIndices.insert(circleIndices.begin(), i);
    }
    faces.push_back(Face(circleIndices));
    points.push_back(top);
    figure.setColor(color);
    figure.setFaces(faces);
    figure.setPoints(points);
    return figure;
}

Figure3D WireFrameParser::parseCylinder(const ini::Configuration &configuration, std::string &name, img::Color &color) {
    double height = configuration[name]["height"].as_double_or_die();
    int n = 2*configuration[name]["n"].as_int_or_die();
    Figure3D figure;
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    std::vector<int> circleIndices = {};

    for(int i = 0; i < n; i ++){
        if (i < n/2) {
            points.push_back(Vector3D::point(cos(4 * i * M_PI / n), sin(4 * i * M_PI / n), 0));
            faces.push_back(Face({i, i + n/2, (i+1)%(n/2) + n/2 , (i+1)%(n/2)}));
            circleIndices.insert(circleIndices.begin(), i);
        }
        if (i == n/2){
            faces.push_back(Face(circleIndices));
            circleIndices = {};
        };
        if(i >= n/2){
            points.push_back(Vector3D::point(cos(4 * i * M_PI / n), sin(4 * i * M_PI / n), height));
            circleIndices.insert(circleIndices.begin(), i);
        }

    }
    faces.push_back(Face(circleIndices));
    figure.setColor(color);
    figure.setPoints(points);
    figure.setFaces(faces);
//    std::cerr << "S: " << circleIndices.size() << std::endl;
    return figure;
}

Figure3D WireFrameParser::parseSphere(const ini::Configuration &configuration, std::string &name, img::Color &color) {
    int n = configuration[name]["n"].as_int_or_die();
    Figure3D figure = parseIcosahedron(color);
    for (int j = 0; j < n; j++) { //repeat the process n times
        std::vector<Vector3D> test = {};
        std::vector<Face> facesTest = {};
        for (Face face1 : figure.getFaces()) {
            Vector3D punt1 = figure.getPoints()[face1.getPointIndices()[0]];
            Vector3D punt2 = figure.getPoints()[face1.getPointIndices()[1]];
            Vector3D punt3 = figure.getPoints()[face1.getPointIndices()[2]];
            Vector3D punt4 = (punt1 + punt2) / 2;
            Vector3D punt5 = (punt2 + punt3) / 2;
            Vector3D punt6 = (punt1 + punt3) / 2;
            int index1, index2, index3, index4, index5, index6;
            if (std::find(test.begin(), test.end(), punt1) == test.end()) {
                index1 = test.size();
                test.push_back(punt1);
            } else {
                index1 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt1));
            }
            if (std::find(test.begin(), test.end(), punt2) == test.end()) {
                index2 = test.size();
                test.push_back(punt2);
            } else {
                index2 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt2));
            }
            if (std::find(test.begin(), test.end(), punt3) == test.end()) {
                index3 = test.size();
                test.push_back(punt3);
            } else {
                index3 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt3));
            }
            if (std::find(test.begin(), test.end(), punt4) == test.end()) {
                index4 = test.size();
                test.push_back(punt4);
            } else {
                index4 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt4));
            }
            if (std::find(test.begin(), test.end(), punt5) == test.end()) {
                index5 = test.size();
                test.push_back(punt5);
            } else {
                index5 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt5));
            }
            if (std::find(test.begin(), test.end(), punt6) == test.end()) {
                index6 = test.size();
                test.push_back(punt6);
            } else {
                index6 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt6));
            }
            facesTest.push_back(Face({index1, index4, index6}));
            facesTest.push_back(Face({index2, index4, index5}));
            facesTest.push_back(Face({index3, index5, index6}));
            facesTest.push_back(Face({index4, index5, index6}));
        }
        figure.setPoints(test);
        figure.setFaces(facesTest);
    }
    std::vector<Vector3D> points = figure.getPoints();
    for (Vector3D &point : points) {
        point /= point.length();
        point.normalise();
    }
    figure.setPoints(points);
    return figure;
}

Figure3D WireFrameParser::parseTorus(const ini::Configuration &configuration, std::string &name, img::Color &color) {
    Figure3D figure;
    double r = configuration[name]["r"].as_double_or_die();
    double R = configuration[name]["R"].as_double_or_die();
    int m = configuration[name]["m"].as_int_or_die();
    int n = configuration[name]["n"].as_int_or_die();
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    for(int i = 0 ; i < n; i++){
        for(int j = 0; j < m; j++){
            double u = 2*i*M_PI/n;
            double v = 2*j*M_PI/m;
            Vector3D point = Vector3D::point((R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v));
            points.push_back(point);
            faces.push_back(Face({m*i+j, m*((i+1)%n) + j, m*((i+1)%n) + (j+1)%m, m*i + (j+1)%m}));
        }
    }
    figure.setColor(color);
    figure.setFaces(faces);
    figure.setPoints(points);
    return figure;
}

Figure3D
WireFrameParser::parse3DLsystem(const ini::Configuration &configuration, std::string &name, img::Color &color) {
    std::string filename = configuration[name]["inputfile"].as_string_or_die();
    LParser::LSystem3D Lsystem = getLSystem3D(filename);
    Figure3D figure = drawLSystem3D(Lsystem, color);
    return figure;
}
//
//Figure3D WireFrameParser::parseMobius(const ini::Configuration &configuration, std::string &name, img::Color &color) {
//    Figure3D figure;
//    double r = configuration[name]["r"].as_double_or_die();
//    double R = configuration[name]["R"].as_double_or_die();
//    int m = configuration[name]["m"].as_int_or_die();
//    int n = configuration[name]["n"].as_int_or_die();
//    std::vector<Vector3D> points;
//    std::vector<Face> faces;
//    for(int i = 0 ; i < n; i++){
//        for(int j = 0; j < m; j++){
//            double u = 2*i*M_PI/n;
//            double v = 2*j*M_PI/m;
//            Vector3D point = Vector3D::point((1+v/2*cos(u/2))*cos(u), (1+v/2*cos(u/2))*sin(u), v/2*sin(u/2));
//            points.push_back(point);
//            faces.push_back(Face({m*i+j, m*((i+1)%n) + j, m*((i+1)%n) + (j+1)%m, m*i + (j+1)%m}));
//        }
//    }
//    figure.setColor(color);
//    figure.setFaces(faces);
//    figure.setPoints(points);
//    return figure;
//}
//
//Figure3D WireFrameParser::parseNavelTorus(const ini::Configuration &configuration, std::string &name, img::Color &color) {
//    Figure3D figure;
//    double r = configuration[name]["r"].as_double_or_die();
//    double R = configuration[name]["R"].as_double_or_die();
//    int m = configuration[name]["m"].as_int_or_die();
//    int n = configuration[name]["n"].as_int_or_die();
//    std::vector<Vector3D> points;
//    std::vector<Face> faces;
//    for(int i = 0 ; i < n; i++){
//        for(int j = 0; j < m; j++){
//            double u = 2*i*M_PI/n;
//            double v = 2*j*M_PI/m;
//            Vector3D point = Vector3D::point(sin(u)*(7+cos(u/3-2*v) + 2* cos(u/3+v)), cos(u)*(7+cos(u/3-2*v) + 2* cos(u/3+v)), sin(u/3-2*v)+2*sin(u/3+v));
//            points.push_back(point);
//            faces.push_back(Face({m*i+j, m*((i+1)%n) + j, m*((i+1)%n) + (j+1)%m, m*i + (j+1)%m}));
//        }
//    }
//    figure.setColor(color);
//    figure.setFaces(faces);
//    figure.setPoints(points);
//    return figure;
//}