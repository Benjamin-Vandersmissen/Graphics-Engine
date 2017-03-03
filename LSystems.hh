//
// Created by Benjamin on 21/02/2017.
//

#ifndef ENGINE_LSYSTEMS_H
#define ENGINE_LSYSTEMS_H
#include "l_parser.hh"
#include "easy_image.hh"
#include <fstream>
#include "Line2D.hh"
#include "UsefulFunctions.hh"

LParser::LSystem2D getLSystem2D(std::string filename);
img::EasyImage
drawLSystem2D(LParser::LSystem2D &lSystem2D, img::Color &bgColor, img::Color &color, int size, bool rainbow);
void LSystem2Dstep(LParser::LSystem2D &lSystem2D, img::Color &color, double &x, double &y, double &angle,
                   Lines2D &lines,
                   int iterations, std::string s, std::vector<std::vector<double>> &brackets);


#endif //ENGINE_LSYSTEMS_H
