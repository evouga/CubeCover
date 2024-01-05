#ifndef TRACESTREAMLINES_H
#define TRACESTREAMLINES_H

#include <Eigen/Core>
#include <vector>

#include "FrameField.h"

namespace CubeCover
{
    class TetMeshConnectivity;
};


struct Streamline
{
    std::vector<int> tetIds;
    std::vector<int> tetFaceIds;
    std::vector<Eigen::Vector3d> directionVecs;
    std::vector<Eigen::Vector3d> points;
    
    std::vector<int> faceIds;

};


void traceStreamlines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const CubeCover::FrameField& frames,
    Eigen::VectorXi& init_tet_ids,
    int max_iter_per_trace,
    std::vector<Streamline>& traces);

#endif