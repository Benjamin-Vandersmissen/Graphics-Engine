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
#include "Figure3D.hh"
#include <algorithm>


LParser::LSystem2D getLSystem2D(std::string filename);
img::EasyImage
drawLSystem2D(LParser::LSystem2D &lSystem2D, img::Color &bgColor, img::Color &color, int size, bool rainbow);
void LSystem2Dstep(LParser::LSystem2D &lSystem2D, img::Color &color, double &x, double &y, double &angle,
                   Lines2D &lines,
                   int iterations, std::string s, std::vector<std::vector<double>> &brackets);

LParser::LSystem3D getLSystem3D(std::string filename);
#endif //ENGINE_LSYSTEMS_H
Figure3D drawLSystem3D(LParser::LSystem3D &lSystem3D, img::Color& color);
void LSystem3Dstep(LParser::LSystem3D &lsystem, std::vector<Face> &faces, Vector3D &point, Vector3D &H, Vector3D &L,
                   Vector3D &U, unsigned int iterations, std::string s, std::vector<std::vector<Vector3D>> &brackets,
                   std::vector<Vector3D> &points);