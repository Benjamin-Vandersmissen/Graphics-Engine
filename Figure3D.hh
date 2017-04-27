//
// Created by uauser on 3/1/17.
//

#ifndef GRAPHICS_ENGINE_FIGURE3D_HH
#define GRAPHICS_ENGINE_FIGURE3D_HH

#include <vector> //std vector
#include "vector.hh" //custom vector
#include "easy_image.hh"
#include "Line2D.hh"
#include <cmath>

class Face {
private:
    /**
     * \brief A vector containing the indices corresponding with the points this face consists of.
     * **/
    std::vector<int> pointIndices;
public:
    /**
     * \brief Constructs a new empty Face.
     * **/
    Face();

    /**
     * \brief Construct a Face from a given vector of indices.
     *
     * \param pointIndices The indices to construct the Face from.
     * **/
    Face(const std::vector<int> &pointIndices);

    /**
     * \brief Returns the pointIndices.
     *
     * \return pointIndices
     * **/
    const std::vector<int> &getPointIndices() const;

    /**
     * \brief Returns the index at position i.
     *
     * \param i The position in pointIndices.
     *
     * \return The index at position i.
     * **/
    unsigned int operator[](unsigned int i) const;
};
/**
 * \brief Triangulates a single Face.
 *
 * \param face The face to triangulate.
 *
 * \return A vector of triangles.
 * **/
std::vector<Face> triangulate(const Face& face);

class Color{
public:
    /***
     * \brief Red color component.
     * */
    double R = 0;
    /**
     * \brief Green color component.
     * **/
    double G = 0;
    /**
     * \brief Blue color component.
     * **/
    double B = 0;

    /**
     * \brief Construct a new Color based on given color components.
     *
     * \param R Red color component.
     *
     * \param G Green color component.
     *
     * \param B Blue color component.
     * **/
    Color(double R, double G, double B);

    /**
     * \brief Constructs a new empty Color.
     * **/
    Color();

    /**
     * \brief Convert to a corresponding img::Color.
     *
     * \return Corresponding img::Color.
     * **/
    img::Color toRGB();

    /**
     * \brief Construct a new Color based on a vector of color components.
     *
     * \param vec An std::vector containing at least three doubles.
     * **/
    Color(std::vector<double> vec);
};

Color operator*(const Color &col1, const Color &col2);

void operator+=(Color& col1, const Color& col2);

Color operator+(const Color& col1, const Color& col2);

Color operator*(const Color& col1, const double d);

class Figure3D {
private:

    /**
     * \brief The faces who form the image.
     * **/
    std::vector<Face> Faces;

    /**
     * \brief The poinst tha make up the image.
     * **/
    std::vector<Vector3D> Points;

    /**
     * \brief The color of the image.
     * **/
    img::Color color;

    /**
     * \brief The reflection rate for ambient light.
     * **/
    Color ambientReflection;

    /**
     * \brief The reflection rate for diffuse light.
     * **/
    Color diffuseReflection;

    /**
     * \brief The reflection rate for specular light.
     * **/
     Color specularReflection;

    /**
     * \brief The reflectioncoefficient.
     * **/
     double reflectionCoefficient;
public:

    /**
     * \brief If the image needs to be drawn in rainbow coloring instead of it's normal color.
     * **/
    bool rainbow = false;

    /**
     * \brief Returns the faces.
     *
     * \return Faces
     * **/
    const std::vector<Face> &getFaces() const;

    /**
     * \brief Sets the faces.
     *
     * \param Faces The new faces.
     * **/
    void setFaces(const std::vector<Face> &Faces);

    /**
     * \brief Returns the points.
     *
     * \return Points
     * **/
    const std::vector<Vector3D> &getPoints() const;

    /**
     * \brief Sets the points.
     *
     * \param Points The new points.
     * **/
    void setPoints(const std::vector<Vector3D> &Points);

    /**
     * \brief Returns the color.
     *
     * \return Color
     * **/
    const img::Color &getColor() const;

    /**
     * \brief Sets the color.
     *
     * \param color The new color.
     * **/
    void setColor(const img::Color &color);

    /**
     * \brief Constructs a new Figure with given faces, points and color.
     *
     * \param color Color
     *
     * \param Points Points
     *
     * \param Faces faces
     * **/
    Figure3D(const std::vector<Face> &Faces, const std::vector<Vector3D> &Points, const img::Color &color);

    /**
     * \brief Constructs a new empty Figure.
     * **/
    Figure3D();

    /**
     * \brief Applies a transformationmatrix to the figure.
     *
     * \param matrix The transformationmatrix.
     * **/
    void applyTransformation(Matrix& matrix);

    /**
     * \brief Compute the center of a given face.
     *
     * \param face The index of the face.
     *
     * \return The centerpoint.
     * **/
    Vector3D getCenter(int face);

    /**
     * \brief Compute the center of the figure.
     *
     * \return The center.
     * **/
    Vector3D getCenter();

    /**
     * \brief Returns a point at a given index.
     *
     * \param i The index of the point.
     *
     * \return The point at index i.
     * **/
    const Vector3D & operator[] (unsigned int i) const;
    /**
     * \brief Returns the ambient reflection component.
     *
     * \return Ambient reflection component.
     * **/
    const Color &getAmbientReflection() const;

    /**
     * \brief Sets the ambient reflection component.
     *
     * \param ambientReflection New ambientReflection.
     * **/
    void setAmbientReflection(const Color &ambientReflection);

    /**
     * \brief Returns the diffuse reflection component.
     *
     * \return Diffuse reflection component.
     * **/
    const Color &getDiffuseReflection() const;

    /**
     * \brief Sets the diffuse reflection component.
     *
     * \param diffuseReflection New diffuseReflection.
     * **/
    void setDiffuseReflection(const Color &diffuseReflection);

    /**
     * \brief Returns the Specular reflection component.
     *
     * \return Specular reflection component.
     * **/
    const Color &getSpecularReflection() const;

    /**
     * \brief Sets the specular reflection component.
     *
     * \param specularReflection New specularReflection.
     * **/
    void setSpecularReflection(const Color &specularReflection);

    /**
     * \brief Returns the reflection coefficient.
     *
     * \return Reflection coefiicient.
     * **/
    double getReflectionCoefficient() const;

    /**
     * \brief Sets the reflection coefficient.
     *
     * \param reflectionCoefficient New reflectionCoefficient.
     * **/
    void setReflectionCoefficient(double reflectionCoefficient);
};

/**
 * \brief Triangulates a figure.
 *
 * \param figure The figure to triangulate.
 * **/
void triangulate(Figure3D& figure);

typedef std::vector<Figure3D> Figures3D;

/**
 * \brief Applies a transformationmatrix to a number of figures.
 *
 * \param figures The figures that will be transformed.
 * **/
void applyTransformation(Figures3D& figures, Matrix& matrix);

/**
 * \brief Constructs a transformationmatrix, based on a given scale.
 *
 * \param scale The scale
 *
 * \return A transformationmatrix
 * **/
Matrix scaleFigure(const double scale);

/**
 * \brief Constructs a transformationmatrix, based on a given rotation.
 *
 * \param angle The angle of the rotation.
 *
 * \return A transformationmatrix
 * */
Matrix rotateFigureX(const double angle);

/**
 * \brief Constructs a transformationmatrix, based on a given rotation.
 *
 * \param angle The angle of the rotation.
 *
 * \return A transformationmatrix
 * */
Matrix rotateFigureY(const double angle);

/**
 * \brief Constructs a transformationmatrix, based on a given rotation.
 *
 * \param angle The angle of the rotation.
 *
 * \return A transformationmatrix
 * */
Matrix rotateFigureZ(const double angle);

/**
 * \brief Constructs a transformationmatrix, based on a given translation.
 *
 * \param vector The tranlation.
 *
 * \return A transformationmatrix
 * **/
Matrix translateFigure(const Vector3D& vector);

/**
 * \brief Constructs a transformation matrix, based on a given point.
 *
 * \param eye The point.
 *
 * \return A transforamtionmatrix
 * **/
Matrix getEyeMatrix(Vector3D& eye);

/**
 * \brief Writes a figure to an output stream.
 *
 * \param stream The outputstream
 *
 * \param figure The figure
 *
 * \return The outputstream
 * **/
std::ostream& operator<<(std::ostream& stream, const Figure3D& figure);

/**
 * \brief Projects 3D figures in a two-dimensional plane.
 *
 * \param figures The figures to be projected.
 *
 * \param ZBuffering Wether or not ZBuffering is applied.
 *
 * \return The projected figures.
 * **/
Lines2D doProjection(const Figures3D &figures, bool ZBuffering);

/**
 * \brief Projects a single point in a two-dimensional plane.
 *
 * \param point The point to be projected.
 *
 * \param d The distance from the eyepoint to the point.
 *
 * \return The projected point.
 * **/
Point2D doProjection(const Vector3D& point, const double d = 1);


/**
 * \brief Compares two faces.
 *
 * \param face1 The first face.
 *
 * \param face2 The second face.
 *
 * \return success
 * **/
bool operator==(const Face& face1, const Face& face2);

#endif //GRAPHICS_ENGINE_FIGURE3D_HH
