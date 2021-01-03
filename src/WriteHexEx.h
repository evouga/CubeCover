#ifndef WRITEHEXEX_H
#define WRITEHEXEX_H

#include <string>
#include <Eigen/Core>

namespace CubeCover
{
    bool writeHexEx(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const Eigen::MatrixXd &vals);
};

#endif