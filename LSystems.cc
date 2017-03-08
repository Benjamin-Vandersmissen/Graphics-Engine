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
    if (iterations > 0) {
        iterations--;
        for (char c : s) {

            if (c != '+' && c != '-' && c != ')' && c != '(') {
                LSystem2Dstep(lSystem2D, color, x, y, angle, lines, iterations,
                              lSystem2D.get_replacement(c), brackets);
            }
            else if (c == '+'){
                angle += toRadial(lSystem2D.get_angle());
            }
            else if (c == '-'){
                angle -= toRadial(lSystem2D.get_angle());
            }
            else if (c == '('){
                brackets.push_back({x,y,angle});
            }else if (c == ')'){
                x = brackets.back()[0];
                y = brackets.back()[1];
                angle = brackets.back()[2];
                brackets.pop_back();
            }
        }
    }else{
        for (char c : s){
            if (c == '+'){
                angle += toRadial(lSystem2D.get_angle());
            }
            else if (c == '-') {
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
            else{
                if (lSystem2D.draw(c)){
                    lines.push_back(Line2D(x,y,x+std::cos(angle), y + std::sin(angle), color));
                }
                x += std::cos(angle);
                y += std::sin(angle);
            }
            if (angle >= toRadial(360)){
                angle -= toRadial(360);
            }
            else if (angle < 0){
                angle += toRadial(360);
            }
        }
        return;
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
    std::vector<Vector3D> points = {};
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

    if (iterations > 0){
        iterations--;
        for(char c : s){
            std::string specialChars = "+-^&\/|()";
            if (specialChars.find(c) == std::string::npos){
                LSystem3Dstep(lsystem, faces, point, H, L, U, iterations,
                              s, brackets, points);
            }
            else if (c == '+'){
                double delta = toRadial(lsystem.get_angle());
                H = H * cos(delta) + L * sin(delta);
                L = -H* sin(delta) + L * cos(delta);
            }
            else if (c == '-'){
                double delta = -toRadial(lsystem.get_angle());
                H = H * cos(delta) + L * sin(delta);
                L = -H* sin(delta) + L * cos(delta);
            }
            else if (c == '^'){
                double delta = toRadial(lsystem.get_angle());
                H = H * cos(delta) + U * sin(delta);
                U = -H*sin(delta) + U * cos(delta);
            }
            else if (c == '&'){
                double delta = -toRadial(lsystem.get_angle());
                H = H * cos(delta) + U * sin(delta);
                U = -H*sin(delta) + U * cos(delta);
            }
            else if (c == '\\'){
                double delta = toRadial(lsystem.get_angle());
                L = L*cos(delta) - U *sin(delta);
                U = L*sin(delta) + U*cos(delta);
            }
            else if (c == '/'){
                double delta = -toRadial(lsystem.get_angle());
                L = L*cos(delta) - U *sin(delta);
                U = L*sin(delta) + U*cos(delta);
            }
            else if (c == '|'){
                H = -H;
                L = -L;
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
            }
        }
    }
    else{
        for(char c : s){
            if (c == '+'){
                double delta = toRadial(lsystem.get_angle());
                H = H * cos(delta) + L * sin(delta);
                L = -H* sin(delta) + L * cos(delta);
            }
            else if (c == '-'){
                double delta = -toRadial(lsystem.get_angle());
                H = H * cos(delta) + L * sin(delta);
                L = -H* sin(delta) + L * cos(delta);
            }
            else if (c == '^'){
                double delta = toRadial(lsystem.get_angle());
                H = H * cos(delta) + U * sin(delta);
                U = -H*sin(delta) + U * cos(delta);
            }
            else if (c == '&'){
                double delta = -toRadial(lsystem.get_angle());
                H = H * cos(delta) + U * sin(delta);
                U = -H*sin(delta) + U * cos(delta);
            }
            else if (c == '\\'){
                double delta = toRadial(lsystem.get_angle());
                L = L*cos(delta) - U *sin(delta);
                U = L*sin(delta) + U*cos(delta);
            }
            else if (c == '/'){
                double delta = -toRadial(lsystem.get_angle());
                L = L*cos(delta) - U *sin(delta);
                U = L*sin(delta) + U*cos(delta);
            }
            else if (c == '|'){
                H = -H;
                L = -L;
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
            }
            else{
                if (lsystem.draw(c)){
                    points.push_back(point);
                }
                point += H;
                if (lsystem.draw(c)){
                    if (points.size() == 1)
                        points.push_back(point);
                    std::vector<int> vector= {points.size()-1, points.size()};
                    faces.push_back(Face(vector));
                }

            }
        }
        return;
    }

}



