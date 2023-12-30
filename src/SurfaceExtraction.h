#ifndef SURFACEEXTRACTION_H
#define SURFACEEXTRACTION_H

#include <Eigen/Core>

namespace CubeCover {

    class TetMeshConnectivity;

    void isosurfaceSoup(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, Eigen::MatrixXd& phivals, Eigen::MatrixXd& isoV, Eigen::MatrixXi& isoF);

    void isosurfaceSoupForSingleIsoVal(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, Eigen::MatrixXd& phival, double isoval, Eigen::MatrixXd& isoV, Eigen::MatrixXi& isoF);

};

#endif