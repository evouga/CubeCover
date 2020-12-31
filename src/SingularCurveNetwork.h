#ifndef SINGULARCURVENETWORK_H
#define SINGULARCURVENETWORK_H

#include <Eigen/Core>

class FrameField;
class TetMeshConnectivity;

namespace CubeCover
{
    void extractSingularCurveNetwork(const Eigen::MatrixXd &V, const TetMeshConnectivity &mesh, const FrameField& field, Eigen::MatrixXd& P, Eigen::MatrixXi& E);
}

#endif