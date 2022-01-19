#ifndef CURLCORRECTION_H
#define CURLCORRECTION_H

#include <Eigen/Core>

class FrameField;

namespace CubeCover
{
    void curlCorrect(const Eigen::MatrixXd& V, FrameField& field, double maxCorrection);
};

#endif