#ifndef EXTRACTISOLINES_H
#define EXTRACTISOLINES_H

#include <Eigen/Core>
#include <vector>
#include <glm/glm.hpp>


namespace CubeCover
{
    class TetMeshConnectivity;
};

void extractIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
    Eigen::MatrixXd& P,
    Eigen::MatrixXi& E,
    Eigen::MatrixXd& P2,
    Eigen::MatrixXi& E2);


void extractPoints(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
    std::vector<glm::vec3>& points,
    std::vector<std::array<double, 3>>& colors );

#endif
