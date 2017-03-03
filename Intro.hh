//
// Created by Benjamin on 14/02/2017.
//

#include "easy_image.hh"
#include "ini_configuration.hh"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef ENGINE_INTRO_H
#define ENGINE_INTRO_H


img::EasyImage ColorRectangle(unsigned int w, unsigned int h, bool scale);
img::EasyImage Blocks(unsigned int Wi, unsigned int Hi, unsigned int nrXBlocks, unsigned int nrYBlocks, std::vector<int> ColorW, std::vector<int> ColorB, bool invert);
img::EasyImage Lines(unsigned int w, unsigned int h, std::string figure, std::vector<int> ColorBG, std::vector<int> ColorLine, unsigned int nrLines);


#endif //ENGINE_INTRO_H
