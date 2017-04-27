//
// Created by Benjamin on 22/02/2017.
//
#include "UsefulFunctions.hh"
int roundToInt(double d){
    return d > 0 ? std::ceil(d-0.5) : std::floor(d + 0.5);

}
std::vector<int> extractColor(std::vector<double> c){
    std::vector<int> color = {(int)(c[0]*255), (int)(c[1]*255), (int)(c[2]*255)};
//    std::cout << color[0] << ", " << color[1] << ", " << color[2] << std::endl;
    color[0] = std::min(255,color[0]);
    color[1] = std::min(255,color[1]);
    color[2] = std::min(255,color[2]);
    return color;
}

double toRadial(double angleInDegrees){
    double angleInRadial = (angleInDegrees/ 360) * (2*3.1415926535897); //approximation
    return angleInRadial;
}

void toPolar(Vector3D &point, double r, double theta, double phi) {
    point.x = r*std::sin(phi)*std::cos(theta);
    point.y = r*std::sin(phi)*std::sin(theta);
    point.z = r*std::cos(phi);
}
