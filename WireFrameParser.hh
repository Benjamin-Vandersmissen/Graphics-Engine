//
// Created by Benjamin on 01/03/2017.
//

#ifndef GRAPHICS_ENGINE_WIREFRAMEPARSER_H
#define GRAPHICS_ENGINE_WIREFRAMEPARSER_H

#include "easy_image.hh"
#include "ini_configuration.hh"
#include "Figure3D.hh"
#include "UsefulFunctions.hh"

class WireFrameParser {
private:
    std::vector<Figure3D> figures;
    img::EasyImage image;
public:
    const img::EasyImage &getImage() const;

public:
    WireFrameParser(const ini::Configuration &configuration);
    Figure3D parseLinedrawing(const ini::Configuration &configuration, std::string &name, img::Color &color);
    Figure3D parseCube(const ini::Configuration &configuration, img::Color &color);
};


#endif //GRAPHICS_ENGINE_WIREFRAMEPARSER_H
