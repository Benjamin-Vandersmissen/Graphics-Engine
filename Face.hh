//
// Created by uauser on 3/1/17.
//

#ifndef GRAPHICS_ENGINE_FACE_HH
#define GRAPHICS_ENGINE_FACE_HH

#include <vector>

class Face {
private:
    std::vector<int> pointIndices;
public:
    Face();

    Face(const std::vector<int> &pointIndices);

    const std::vector<int> &getPointIndices() const;
};


#endif //GRAPHICS_ENGINE_FACE_HH
