//
// Created by Benjamin on 01/03/2017.
//

#ifndef GRAPHICS_ENGINE_WIREFRAMEPARSER_H
#define GRAPHICS_ENGINE_WIREFRAMEPARSER_H

#include "easy_image.hh"
#include "ini_configuration.hh"
#include "Figure3D.hh"
#include "UsefulFunctions.hh"
#include "LSystems.hh"
#include <algorithm>

class WireFrameParser {
private:
    std::vector<Figure3D> figures;
    img::EasyImage image;
public:
    const img::EasyImage &getImage() const;

public:
    WireFrameParser(const ini::Configuration &configuration);
    Figure3D parseLinedrawing(const ini::Configuration &configuration, std::string &name, img::Color &color);
    Figure3D parseCube(img::Color &color);
    Figure3D parseTetrahedron(img::Color& color);
    Figure3D parseOctahedron(img::Color& color);
    Figure3D parseIcosahedron(img::Color& color);
    Figure3D parseDodecahedron(img::Color& color);
    Figure3D parseCone(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseCylinder(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseSphere(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseTorus(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parse3DLsystem(const ini::Configuration &configuration, std::string & name, img::Color& color);
};


#endif //GRAPHICS_ENGINE_WIREFRAMEPARSER_H
