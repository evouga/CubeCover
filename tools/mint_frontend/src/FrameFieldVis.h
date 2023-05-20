#ifndef FRAMEFIELDVIS_H
#define FRAMEFIELDVIS_H

#include <Eigen/Core>
#include <vector>

namespace CubeCover
{
    class TetMeshConnectivity;
    class FrameField;
};

void buildFrameVectors(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    double scale,
    Eigen::MatrixXd& centroids,
    std::vector<Eigen::MatrixXd>& frameVectors
    );


void computePerVectorCurl(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    Eigen::MatrixXd& splitCurl
    );    

#endif