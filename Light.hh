//
// Created by Benjamin on 25/04/2017.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "Figure3D.hh"
#include "vector.hh"

class Light {
public:
    /**
     * \brief The ambient light component.
     * **/
    Color ambientLight;
    /**
     * \brief The diffuse light component.
     * **/
    Color diffuseLight;
    /**
     * \brief The specular light component.
     * **/
    Color specularLight;

    /**
     * \brief Infinite lightsource or point?
     * **/
    bool infinity = true;

    /**
     * \brief The location of the lightsource when it's not infinite.
     * **/
    Vector3D location = Vector3D::point(0,0,0);

    /**
     * \brief The direction of the lightsource when it's infinite.
     * **/
    Vector3D direction = Vector3D::vector(0,0,0);
};


typedef std::vector<Light> Lights3D;

#endif //ENGINE_LIGHT_H
