#pragma once

#include <Eigen/Core>
#include "TetMeshConnectivity.h"

namespace CubeCover {
  // Extract isolines of the given values (ntetx x 3)
  void ExtractIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
                       Eigen::MatrixXd& P,
                       Eigen::MatrixXi& E,
                       Eigen::MatrixXd& P2,
                       Eigen::MatrixXi& E2);
}