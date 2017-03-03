//
// Created by Benjamin on 14/02/2017.
//

#include "easy_image.hh"
#include "ini_configuration.hh"
#include "UsefulFunctions.hh"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef ENGINE_INTRO_H
#define ENGINE_INTRO_H


class IntroParser{
private:
    img::EasyImage image;
protected:
    void parseColorRectangle(const ini::Configuration &configuration);
    void parseBlocks(const ini::Configuration &configuration);
    void parseLines(const ini::Configuration &configuration);
    void QuarterCircle(unsigned int w, unsigned int h, std::vector<int> ColorLine, unsigned int nrLines,
                           unsigned int quadrant = 1, unsigned int x = 0, unsigned int y = 0);
    void Eye(unsigned int w, unsigned int h, std::vector<int> ColorLine, unsigned int nrLines);
    void Diamond(unsigned int w, unsigned int h, std::vector<int> ColorLine, unsigned int nrLines);
public:
    IntroParser(const ini::Configuration &configuration);
    const img::EasyImage &getImage() const;
};



#endif //ENGINE_INTRO_H
