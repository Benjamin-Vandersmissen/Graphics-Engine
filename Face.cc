//
// Created by uauser on 3/1/17.
//

#include "Face.hh"

Face::Face(const std::vector<int> &pointIndices) : pointIndices(pointIndices) {}

Face::Face() {}

const std::vector<int> &Face::getPointIndices() const {
    return pointIndices;
}
