#ifndef VOXELUTILS_H
#define VOXELUTILS_H

#include <Eigen/Core>
#include <vector>
#include <glm/glm.hpp>

bool pointInsideT(const Eigen::Vector3d& A, 
	              const Eigen::Vector3d& B, 
	              const Eigen::Vector3d& C, 
	              const Eigen::Vector3d& D, 
	              const Eigen::Vector3d& samplePoint, Eigen::VectorXd& textureCoordinate);

/*

namespace CubeCover
{
    class TetMeshConnectivity;
};

void extractIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
    Eigen::MatrixXd& P,
    Eigen::MatrixXi& E,
    Eigen::MatrixXd& P2,
    Eigen::MatrixXi& E2,
    std::vector<int> badverts);


void extractPoints(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
    std::vector<glm::vec3>& points,
    std::vector<std::array<double, 3>>& colors );
*/

#endif