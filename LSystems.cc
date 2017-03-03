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
//        std::cout << s << std::endl;
        for (char c : s) {

//            std::cout << "Iteration" << iterations << std::endl;
//            std::cerr << x << ',' << y << std::endl;
            if (c != '+' && c != '-' && c != ')' && c != '(') {
                LSystem2Dstep(lSystem2D, color, x, y, angle, lines, iterations,
                              lSystem2D.get_replacement(c), brackets);
//                std::cout << "~~~" << std::endl;
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
//            std::cout << c << std::endl;
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
//                std::cerr << angle << std::endl;

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
//    for(Line2D line :lines){
//        std::cerr << line.point1.x << ' ' << line.point1.y << ' ' << line.point2.x << ' ' << line.point2.y << std::endl;
//    }
    img::EasyImage image = draw2DLines(lines, size, bgColor, rainbow);
    return image;
}



