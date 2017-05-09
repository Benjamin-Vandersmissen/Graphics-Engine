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
#include "Light.hh"
#include <algorithm>

class WireFrameParser {
private:
    std::vector<Figure3D> figures;
    img::EasyImage image;
    unsigned int ZBuffering;
protected:
    void rearrangeTriangles(std::vector<Face> &triangles, std::vector<Vector3D> points);
    Figure3D drawCube(Color ambientReflection, Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient,
                          Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawTetrahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                 unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawOctahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawIcosahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                 unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D
    drawDodecahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                         unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawCone(double height, int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                          unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawCuboid(double length, double height, double depth, Color ambientReflection, Color diffuseReflection,
                            Color specularReflection, unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0),
                            double scale = 1);
    Figure3D drawCylinder(double height, int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                              unsigned int reflectionCoefficient, bool zijvlakken = true, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawSphere(int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                            unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawTorus(double r, double R, int m, int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                           unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawBuckyBall(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                               unsigned int reflectionCoefficient, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figure3D drawMengerSponge(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                  unsigned int reflectionCoefficient, const int iterations, Vector3D center = Vector3D::point(0,0,0), Vector3D rotation = Vector3D::vector(0,0,0), double scale = 1);
    Figures3D generateFractal(Figure3D& figure, const int iterations, const double scale);

    Figures3D makeThicc(Figure3D &figure, const double radius, const int n, const int m);

public:
    const img::EasyImage &getImage() const;
    WireFrameParser(const ini::Configuration &configuration, unsigned int ZBuffering = 0);
    Figure3D parseLinedrawing(const ini::Configuration &configuration, std::string &name);
    Figure3D parseCube(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                           unsigned int reflectionCoefficient);
    Figure3D parseTetrahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                  unsigned int reflectionCoefficient);
    Figure3D parseOctahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                 unsigned int reflectionCoefficient);
    Figure3D parseIcosahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                  unsigned int reflectionCoefficient);
    Figure3D parseDodecahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                   unsigned int reflectionCoefficient);
    Figure3D parseCone(const ini::Configuration &configuration, std::string &name, Color ambientReflection, Color diffuseReflection,
                           Color specularReflection, unsigned int reflectionCoefficient);
    Figure3D parseCuboid(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                             Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient);
    Figure3D parseCylinder(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                               Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient);
    Figure3D parseSphere(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                             Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient);
    Figure3D parseTorus(const ini::Configuration &configuration, std::string &name, Color ambientReflection, Color diffuseReflection,
                            Color specularReflection, unsigned int reflectionCoefficient);
    Figure3D parse3DLsystem(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient);
    Figures3D parseFractal(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                               Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient);
    Figure3D parseBuckyBall(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                unsigned int reflectionCoefficient);
    Figure3D parseMengerSponge(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                   Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient);
    Figures3D parseThick(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                         Color diffuseReflection,
                         Color specularReflection, unsigned int reflectionCoefficient);

    /**
     * @brief draws a railroad track from input
     * **/

    std::vector<Figure3D> parseRail();
    std::vector<Figure3D> parseTrain(Color baseColor);
    std::vector<Figure3D> parseStation(Color baseColor);
    std::vector<Figure3D> parseDirections();
    Figure3D mergeFigures(Figures3D& figures);
//    Figure3D parseMobius(const ini::Configuration & configuration, std::string& name, img::Color& color);
//    Figure3D parseNavelTorus(const ini::Configuration &configuration, std::string &name, img::Color& color);

};

/**
 * \brief Draw a single ZBuffered triangle.
 *
 * \param buf The ZBuffer
 *
 * \param image The image to draw on.
 *
 * \param A The first point of the triangle.
 *
 * \param B The second point of the triangle.
 *
 * \param C The third point of the triangle.
 *
 * \param d The distance from the eye to the triangle.
 *
 * \param dx The amount by which the projected triangle needs to be moved in the x-plane.
 *
 * \param dy The amount by which the projected triangle needs to be moved in the y-plane.
 *
 * \param color The color
 * **/
void draw_zbuf_triangle(ZBuffer &buf, img::EasyImage &image, Vector3D &A, Vector3D &B, Vector3D &C, double d, double dx,
                        double dy, Lights3D &lights, const Color &ambientReflection, const Color &diffuseReflection,
                        const Color &specularReflection, double reflectionCoefficient, Vector3D eye);

void draw_textured_triangle(ZBuffer &buf, img::EasyImage &image, Vector3D &A, Vector3D &B, Vector3D &C, double d, double dx,
                            double dy, std::vector<Vector3D> rectangleProperties, img::EasyImage *texture, Vector3D eye);

#endif //GRAPHICS_ENGINE_WIREFRAMEPARSER_H
