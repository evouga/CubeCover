#ifndef WRITEFRAMEFIELD_H
#define WRITEFRAMEFIELD_H

#include <Eigen/Core>

namespace CubeCover
{
    class FrameField;

    bool writeFrameField(const std::string& fraFilename, const std::string& permFilename, const FrameField &field);

};

#endif