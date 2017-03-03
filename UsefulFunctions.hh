//
// Created by Benjamin on 22/02/2017.
//

#ifndef ENGINE_USEFULFUNCTIONS_H
#define ENGINE_USEFULFUNCTIONS_H
#include <cmath>
#include <vector>
#include "vector.hh"
int roundToInt(double d);
std::vector<int> extractColor(std::vector<double> c);
double toRadial(double angleInDegrees);
void toPolar(Vector3D &point, double r, double theta, double phi);

#endif //ENGINE_USEFULFUNCTIONS_H
