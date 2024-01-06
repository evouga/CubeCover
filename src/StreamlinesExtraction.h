#pragma once

#include <Eigen/Core>
#include <vector>

#include "FrameField.h"
#include "TetMeshConnectivity.h"

namespace CubeCover {

struct Streamline
{
  std::vector<int> tetIds;
  std::vector<int> tetFaceIds;
  std::vector<Eigen::Vector3d> directionVecs;
  std::vector<Eigen::Vector3d> points;

  std::vector<int> faceIds;

};

// Initialize the tets which may contain one or more isolines
std::vector<int> InitializeTracingTets(const TetMeshConnectivity& mesh, const Eigen::MatrixXd& values);

// Tracing stream lines
void TraceStreamlines(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const FrameField& frames,
                      std::vector<int>& init_tet_ids,
                      int max_iter_per_trace,
                      std::vector<Streamline>& traces);
}