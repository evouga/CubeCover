#pragma once

#include <Eigen/Core>
#include <vector>

#include "TetMeshConnectivity.h"

namespace CubeCover {

// The streamline tracing point: start from the start_pt_ and the tracing direction is given by dir_
// To identify whether this stream pts is redundant: two stream pts has the "same" start_pts and same dir
struct StreamPt {
  StreamPt(const Eigen::Vector3d start_pt, const Eigen::Vector3d& dir, double eps = 1e-6, int face_id = -1) : start_pt_(start_pt), dir_(dir), face_id_(face_id), eps_(eps) {
    dir_.normalize();
  }

  bool operator==(const StreamPt& other) const {
    return (this->start_pt_ - other.start_pt_).norm() <= eps_ && (this->dir_ - other.dir_).norm() == 0 && (this->face_id_ == other.face_id_);
  }

  // a really naive hash function
  struct HashFunction {
    size_t operator()(const StreamPt& stream_pt) const {
      std::size_t h1 = std::hash<double>{}(stream_pt.dir_[0]);
      std::size_t h2 = std::hash<double>{}(stream_pt.dir_[1]);
      std::size_t h3 = std::hash<double>{}(stream_pt.dir_[2]);

      std::size_t h4 = std::hash<int>{}(int(std::floor(stream_pt.start_pt_[0])));
      std::size_t h5 = std::hash<int>{}(int(std::floor(stream_pt.start_pt_[1])));
      std::size_t h6 = std::hash<int>{}(int(std::floor(stream_pt.start_pt_[2])));
      std::size_t h7 = std::hash<int>{}(stream_pt.face_id_);
      return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4) ^ (h6 << 5) ^ (h7 << 6);
    }
  };

  Eigen::Vector3d start_pt_;
  Eigen::Vector3d dir_;
  int face_id_;                     // optional, specify which face this point lies on. We only eliminate the stream pts which are close and on the same face
  double eps_;
};

// The streamline tracing Tet and tracing dierection. This is used to prevent the loop during the tracing
struct StreamTetDir {
    StreamTetDir(int tet_id, const Eigen::Vector3d& dir) : dir_(dir), tet_id_(tet_id) {
        dir_.normalize();
    }

    bool operator==(const StreamTetDir& other) const {
        return (this->tet_id_ == other.tet_id_) && (this->dir_ - other.dir_).norm() == 0;
    }

    struct HashFunction {
        size_t operator()(const StreamTetDir& tet_dir) const {
            std::size_t h1 = std::hash<double>{}(tet_dir.dir_[0]);
            std::size_t h2 = std::hash<double>{}(tet_dir.dir_[1]);
            std::size_t h3 = std::hash<double>{}(tet_dir.dir_[2]);

            std::size_t h4 = std::hash<int>{}(tet_dir.tet_id_);

            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
        }
    };

    Eigen::Vector3d dir_;
    int tet_id_;
};

// The streamline data structure
struct Streamline {
  std::vector<int> tet_ids_;                       // the tet ids which this streamline passes
  std::vector<int> tet_face_ids_;                  // the corresponding tet face ids, where the streamline intersect with the tet
  std::vector<StreamPt> stream_pts_;               // the stream pts. If tet_face_ids[i] != -1, then the stream_pts_.start_pt_ is on that face. Otherwise, it is inside the tet_ids[i]
};

// check whether a line: l(t) = p0 + t * dir is intersection with the triangle v0, v1, v2;
// If intersection happens at p = p0 + t * dir, we have p - v0 = c0 (v1 - v0) + c1 (v2- v0) => c0 * (v1 - v0)  + c1 * (v2 - v0) - t * dir = p0 - v0
// We can solve the linear system [v1 - v0, v2 - v0, -dir] [c0, c1, t]^T = [p0 - v0]. This can be solved if dir is not on the triangle plane
// Remark: if you want to check the intersection between ray and the triangle, you can add the check for t >= 0.
bool IntersectTriangle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& p0, const Eigen::Vector3d& dir, double& t, bool is_debug);

// Random Sampling the stream points: sampling "nsamples" tets randomly as the seed of the streamline tracing. If is_random_inside = false, we also random sample a point inside the sampled the tet.
// Otherwise, we use the tet centroid
std::vector<std::pair<int, Eigen::Vector3d>> RandomSampleStreamPts(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, int nsamples, bool is_random_inside = false);

// Compute the integer grid point within the tet. Assume p = c0 V0 + c1 V1 + c2 V2 + c3 V3. For the j-th channel of values, we have \sum ci * val_ij = k. If we have 3 channels, then we have unique solution.
// Remark: values are defined on the tet soup, that is for tet i, tet_vertex id j, its value is equal to values(4 * i + j) NOT values(mesh.tetVertex(i, j))
std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> ComputeGridPts(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& values, int tet_id);

// Tracing streamlines
void TraceStreamlines(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& frames, int max_iter_per_trace,
                      std::vector<Streamline>& traces, int num_seeds = 100, bool is_random_inside = false, double eps = 1e-6);

// Tracing streamlines from integer points. 
void TraceStreamlines(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& frames, const Eigen::MatrixXd& values, int max_iter_per_trace,
    std::vector<Streamline>& traces, double eps = 1e-6);

// Tracing streamlines if the tracing seeds (stream pts) are given
void TraceStreamlines(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& frames, const std::vector<std::pair<int, StreamPt>>& init_tracing_pts, int max_iter_per_trace,
    std::vector<Streamline>& traces);

// DFS helper
bool AdvanceStreamline(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& frames,
                       Streamline& s, bool is_debug = false);
}