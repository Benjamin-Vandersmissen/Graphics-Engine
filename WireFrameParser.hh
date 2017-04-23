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
    unsigned int ZBuffering;
protected:
    Figure3D drawCube(img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawTetrahedron(img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawOctahedron(img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawIcosahedron(img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawDodecahedron(img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawCone(double height, int n, img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawCuboid(double height, double length, double depth, img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawCylinder(double height, int n, img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawSphere(int n, img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawTorus(double r, double R, int m, int n, img::Color &color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawBuckyBall(img::Color& color, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawMengerSponge(img::Color &color, const int iterations, Vector3D center = Vector3D::point(0, 0, 0),
                              Vector3D rotation = Vector3D::vector(0, 0, 0), double scale = 1);
    Figures3D generateFractal(Figure3D& figure, const int iterations, const double scale);

public:
    const img::EasyImage &getImage() const;
    WireFrameParser(const ini::Configuration &configuration, unsigned int ZBuffering = 0);
    Figure3D parseLinedrawing(const ini::Configuration &configuration, std::string &name, img::Color &color);
    Figure3D parseCube(img::Color &color);
    Figure3D parseTetrahedron(img::Color& color);
    Figure3D parseOctahedron(img::Color& color);
    Figure3D parseIcosahedron(img::Color& color);
    Figure3D parseDodecahedron(img::Color& color);
    Figure3D parseCone(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseCuboid(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseCylinder(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseSphere(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parseTorus(const ini::Configuration &configuration, std::string &name, img::Color& color);
    Figure3D parse3DLsystem(const ini::Configuration &configuration, std::string & name, img::Color& color);
    Figures3D parseFractal(const ini::Configuration& configuration, std::string& name, img::Color& color);
    Figure3D parseBuckyBall(img::Color& color);
    Figure3D parseMengerSponge(const ini::Configuration &configuration, std::string &name, img::Color &color);

    /**
     * @brief draws a railroad track from input
     * **/

    std::vector<Figure3D> parseRail();
    std::vector<Figure3D> parseTrain();
    std::vector<Figure3D> parseDirections();
    Figure3D mergeFigures(Figures3D& figures);
//    Figure3D parseMobius(const ini::Configuration & configuration, std::string& name, img::Color& color);
//    Figure3D parseNavelTorus(const ini::Configuration &configuration, std::string &name, img::Color& color);

};


#endif //GRAPHICS_ENGINE_WIREFRAMEPARSER_H
