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


void makeEdgeSpanningTree(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    int startTetId,
    std::vector<Eigen::Vector2i>& tree_traversal,
    std::vector<Eigen::Vector4i>& tree_traversal_metadata
);

void integrateFieldOnEdges(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    const std::vector<Eigen::Vector2i>& tree_traversal,
    const std::vector<Eigen::Vector4i>& tree_traversal_metadata,
    double period,
    Eigen::MatrixXd& integratedVals
);

void projectVertScalarsToTetFrames(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    const Eigen::MatrixXd& integratedVals,
    std::vector<Eigen::MatrixXd>& ret_frames
);


#endif