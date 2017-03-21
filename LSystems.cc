//
// Created by Benjamin on 21/02/2017.
//

#include "LSystems.hh"
LParser::LSystem2D getLSystem2D(std::string filename) {
    std::ifstream stream(filename);
    LParser::LSystem2D lSystem2D(stream);
    return lSystem2D;
}

void LSystem2Dstep(LParser::LSystem2D &lSystem2D, img::Color &color, double &x, double &y, double &angle,
                   Lines2D &lines,
                   int iterations, std::string s, std::vector<std::vector<double>> &brackets) {
    for (char c : s) {
        std::string specialChars = "+-()";
        if (specialChars.find(c) == std::string::npos) {
            if (iterations > 0) {
                LSystem2Dstep(lSystem2D, color, x, y, angle, lines, iterations - 1,
                              lSystem2D.get_replacement(c), brackets);
            }
            else{
                if (lSystem2D.draw(c)){
                    lines.push_back(Line2D(x,y,x+std::cos(angle), y + std::sin(angle), color));
                }
                x += std::cos(angle);
                y += std::sin(angle);
            }
        }
        else if (c == '+'){
            angle += toRadial(lSystem2D.get_angle());
        }
        else if (c == '-'){
            angle -= toRadial(lSystem2D.get_angle());
        }
        else if (c == '('){
            brackets.push_back({x,y,angle});
        }
        else if (c == ')'){
            x = brackets.back()[0];
            y = brackets.back()[1];
            angle = brackets.back()[2];
            brackets.pop_back();
        }
    }
}

img::EasyImage drawLSystem2D(LParser::LSystem2D &lSystem2D, img::Color &bgColor, img::Color &color, int size, bool rainbow = false) {
    Lines2D lines;
    std::string s = lSystem2D.get_initiator();
    double x = 0;
    double y = 0;
    double angle = lSystem2D.get_starting_angle();
    std::vector<std::vector<double>> brackets;
    LSystem2Dstep(lSystem2D, color, x, y, angle, lines, lSystem2D.get_nr_iterations(),
                  lSystem2D.get_initiator(), brackets);

    img::EasyImage image = draw2DLines(lines, size, bgColor, rainbow);
    return image;
}

LParser::LSystem3D getLSystem3D(std::string filename) {
    std::ifstream file(filename);
    LParser::LSystem3D lsystem(file);
    return lsystem;
}

Figure3D drawLSystem3D(LParser::LSystem3D &lSystem3D, img::Color &color) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(0,0,0)};
    std::vector<Face> faces = {};
    std::string s = lSystem3D.get_initiator();
    Vector3D point = Vector3D::point(0,0,0);
    Vector3D H = Vector3D::vector(1,0,0);
    Vector3D L = Vector3D::vector(0,1,0);
    Vector3D U = Vector3D::vector(0,0,1);
    std::vector<std::vector<Vector3D>> brackets = {};

    LSystem3Dstep(lSystem3D, faces, point, H, L, U, lSystem3D.get_nr_iterations(),
                  s, brackets, points);


    figure.setColor(color);
    figure.setFaces(faces);
    figure.setPoints(points);
    return figure;
}

void LSystem3Dstep(LParser::LSystem3D &lsystem, std::vector<Face> &faces, Vector3D &point, Vector3D &H, Vector3D &L,
                   Vector3D &U, unsigned int iterations, std::string s, std::vector<std::vector<Vector3D>> &brackets,
                   std::vector<Vector3D> &points) {
    for(char c : s){
//            std::cerr << "H: " << H << std::endl;
//            std::cerr << "L: " << L << std::endl;
//            std::cerr << "U: " << U << std::endl;
        std::string specialChars = "+-^&\\/|()";
        if (specialChars.find(c) == std::string::npos) {
            if (iterations > 0) {
                LSystem3Dstep(lsystem, faces, point, H, L, U, iterations-1,
                              lsystem.get_replacement(c), brackets, points);
            }
            else{
                point += H;
                if (lsystem.draw(c)){
                    points.push_back(point);
                    std::vector<int> indices = {points.size()-2, points.size()-1};
                    faces.push_back(Face({indices}));
                }
            }
        }
        else if (c == '+'){
            double delta = toRadial(lsystem.get_angle());
            Vector3D newH = H*cos(delta) + L*sin(delta);
            Vector3D newL = -H*sin(delta) + L*cos(delta);
            H = newH;
            L = newL;
        }
        else if (c == '-'){
            double delta = -toRadial(lsystem.get_angle());
            Vector3D newH = H*cos(delta) + L*sin(delta);
            Vector3D newL = -H*sin(delta) + L*cos(delta);
            H = newH;
            L = newL;
        }
        else if (c == '^'){
            double delta = toRadial(lsystem.get_angle());
            Vector3D newH = H * cos(delta) + U * sin(delta);
            Vector3D newU = -H*sin(delta) + U * cos(delta);
            H = newH;
            U = newU;
        }
        else if (c == '&'){
            double delta = -toRadial(lsystem.get_angle());
            Vector3D newH = H * cos(delta) + U * sin(delta);
            Vector3D newU = -H*sin(delta) + U * cos(delta);
            H = newH;
            U = newU;
        }
        else if (c == '\\'){
            double delta = toRadial(lsystem.get_angle());
            Vector3D newL = L*cos(delta) - U*sin(delta);
            Vector3D newU = L*sin(delta) + U*cos(delta);
            L = newL;
            U = newU;
        }
        else if (c == '/'){
            double delta = -toRadial(lsystem.get_angle());
            Vector3D newL = L*cos(delta) - U*sin(delta);
            Vector3D newU = L*sin(delta) + U*cos(delta);
            L = newL;
            U = newU;
        }
        else if (c == '|'){
            Vector3D newH = -H;
            Vector3D newL = -L;
            H = newH;
            L = newL;
        }
        else if (c == '('){
            brackets.push_back({point, H, L, U});
        }
        else if (c == ')'){
            point = brackets.back()[0];
            H = brackets.back()[1];
            L = brackets.back()[2];
            U = brackets.back()[3];
            brackets.pop_back();
            points.push_back(point);
        }
    }
}



