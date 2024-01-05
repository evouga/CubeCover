#ifndef EXTRACTISOLINES_H
#define EXTRACTISOLINES_H

#include <Eigen/Core>

namespace CubeCover
{
    class TetMeshConnectivity;
};

void extractIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
                     Eigen::MatrixXd& P,
                     Eigen::MatrixXi& E,
                     Eigen::MatrixXd& P2,
                     Eigen::MatrixXi& E2);

#endif