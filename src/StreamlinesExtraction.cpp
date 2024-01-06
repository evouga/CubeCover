#include "StreamlinesExtraction.h"

#include <cassert>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <random>

namespace CubeCover {
// check whether a line: l(t) = p0 + t * dir is intersection with the triangle v0, v1, v2;
bool IntersectTriangle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& p0, const Eigen::Vector3d& dir, double& t, bool is_debug) {
  Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
  n.normalize();
  
  // parallel case
  if(std::abs(n.dot(dir)) < 1e-6) {
    if(is_debug) {
      std::cout << "parallel case: " << std::endl;
      std::cout << "normal: \n" << n << std::endl;
      std::cout << "dir: \n" << dir << std::endl;
    }
    return false;
  }

  Eigen::Matrix3d A;
  A.col(0) = v1 - v0;
  A.col(1) = v2 - v0;
  A.col(2) = -dir;

  Eigen::Vector3d rhs = p0 - v0;

  if(std::abs(A.determinant()) < 1e-10) {
    return false;
  }
  Eigen::Vector3d sol = A.inverse() * rhs;

  t = sol[2];
  return (sol[0] >= 0) && (sol[1] >= 0) && (sol[0] + sol[1]) <= 1;
}

// return a random number between [min, max]
static double RandomInRange(double min, double max) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  return min + (max - min) * dis(gen);
}

// sample a point inside the tet (may lay on the boundary): if not random, return the centroid
static Eigen::Vector3d SamplePointInsideTet(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, int tet_id, bool is_random = false) {
  std::array<double, 4> barycentrics = {0.25, 0.25, 0.25, 0.25};
  if(is_random) {
    barycentrics[0] = RandomInRange(0, 1);
    barycentrics[1] = RandomInRange(0, 1 - barycentrics[0]);
    barycentrics[2] = RandomInRange(0, 1 - barycentrics[0] - barycentrics[1]);
    barycentrics[3] = 1 - barycentrics[0] - barycentrics[1] - barycentrics[2];
  }

  Eigen::Vector3d pt = Eigen::Vector3d::Zero();
  for(int i = 0; i < 4; i++) {
    pt += V.row(mesh.tetVertex(tet_id, i)).transpose() * barycentrics[i];
  }
  return pt;
}

// Random Sampling the stream points: sampling "nsamples" tets randomly as the seed of the streamline tracing. If is_random_inside = false, we also random sample a point inside the sampled the tet.
// Otherwise, we use the tet centroid
std::vector<std::pair<int, Eigen::Vector3d>> RandomSampleStreamPts(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, int nsamples, bool is_random_inside) {
  int ntets = mesh.nTets();
  std::vector<std::pair<int, Eigen::Vector3d>> samples;

  if(ntets <= nsamples) {
    samples.resize(ntets);
    for(int i = 0; i < ntets; i++) {
      samples[i].first = i;
      samples[i].second = SamplePointInsideTet(V, mesh, i, is_random_inside);
    }
    return samples;
  } else {
    std::unordered_set<int> sampled_ids;
    while(sampled_ids.size() < nsamples) {
      int id = std::rand() % ntets;
      if(!sampled_ids.count(id)) {
        sampled_ids.insert(id);
      }
    }
    std::vector<int> tmp_vec;
    tmp_vec.insert(tmp_vec.end(), sampled_ids.begin(), sampled_ids.end());

    for(auto tid : tmp_vec) {
      Eigen::Vector3d pt = SamplePointInsideTet(V, mesh, tid, is_random_inside);
      samples.push_back({tid, pt});
    }
  }
  return samples;
}

// Check whether there is f(x) = iso_val across this segment whose endpoints' values are val1 and val2.
static bool IsCrossing(double val1, double val2, double iso_val, bool flip) {
  if(!flip) {
    return val1 < iso_val && iso_val <= val2;
  } else {
    return val1 <= iso_val && iso_val < val2;
  }
}

// Initialize the tets which may contain one or more isolines
std::vector<int> InitializeTracingTets(const TetMeshConnectivity& mesh, const Eigen::MatrixXd& values) {
  int ntets = mesh.nTets();

  auto has_cross = [&mesh, &values](int tet_id, int channel) {
    double minphi = std::numeric_limits<double>::infinity();
    double maxphi = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < 4; i++){
      double phival = values(4 * tet_id + i, channel);   // value is defined on the tet soup!
      minphi = std::min(phival, minphi);
      maxphi = std::max(phival, maxphi);
    }

    int minint = int(minphi);
    int maxint = int(maxphi);
    bool is_cross = false;
    for (int isoval = minint; isoval <= maxint && !is_cross; isoval++) {
      for (int idx1 = 0; idx1 < 4 && !is_cross; idx1++) {
        for (int idx2 = idx1 + 1; idx2 < 4 && !is_cross; idx2++) {
          int vidx1 = mesh.tetVertex(tet_id, idx1);
          int vidx2 = mesh.tetVertex(tet_id, idx2);
          double val1 = values(4 * tet_id + idx1, channel);
          double val2 = values(4 * tet_id + idx2, channel);
          if (val1 > val2) {
            std::swap(val1, val2);
            std::swap(vidx1, vidx2);
          }
          if (IsCrossing(val1, val2, double(isoval), 0) ||
              IsCrossing(val1, val2, double(isoval), 1)) {
            is_cross = true;
          }
        }
      }
    }
    return is_cross;
  };
  std::vector<int> tet_ids;
  for(int i = 0; i < ntets; i++) {
    bool is_cross0 = has_cross(i, 0);
    bool is_cross1 = has_cross(i, 1);
    bool is_cross2 = has_cross(i, 2);
    if((is_cross0 && is_cross1) || (is_cross0 && is_cross2) || (is_cross1 && is_cross2)) {
      tet_ids.push_back(i);
    }
  }
  return tet_ids;
}


bool AdvanceStreamline(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& frames,
                       Streamline& s, double stream_pt_eps, bool is_debug) {
  int cur_tet_id = s.tet_ids_.back();
  int nvecs = frames.rows() / mesh.nTets();

  Eigen::MatrixXd cur_tet_frame = frames.block(cur_tet_id * nvecs, 0, nvecs, 3);
  Eigen::Vector3d cur_point = s.stream_pts_.back().start_pt_;
  Eigen::Vector3d prev_vec = s.stream_pts_.back().dir_;

  double max_match = 0;
  Eigen::Vector3d first_frame_vec = cur_tet_frame.row(0);
  Eigen::Vector3d cur_vec = prev_vec.cross(first_frame_vec);
  for (int i = 0; i < nvecs; i++) {
    Eigen::Vector3d cur_dir = cur_tet_frame.row(i);
    double cur_match = cur_dir.dot(prev_vec) / cur_dir.norm() / prev_vec.norm();

    if (std::abs(cur_match) > max_match) {
      max_match = std::abs(cur_match);
      cur_vec = cur_tet_frame.row(i);
      if (cur_match < 0)
        cur_vec = -cur_vec;
    }
  }
  cur_vec.normalize();

  int cur_local_fid = s.tet_face_ids_.back();

  int next_tetface_id = -1;
  int next_face_id = -1;
  int next_tet_id = -1;
  Eigen::Vector3d intersect_point = cur_point;
  int next_face_orientation = -1;

  for (int i = 0; i < 4; i++) {
    if (i != cur_local_fid) {
      int fid = mesh.tetFace(cur_tet_id, i);
      Eigen::Vector3d v0, v1, v2;
      v0 << V(mesh.faceVertex(fid, 0), 0), V(mesh.faceVertex(fid, 0), 1),
          V(mesh.faceVertex(fid, 0), 2);
      v1 << V(mesh.faceVertex(fid, 1), 0), V(mesh.faceVertex(fid, 1), 1),
          V(mesh.faceVertex(fid, 1), 2);
      v2 << V(mesh.faceVertex(fid, 2), 0), V(mesh.faceVertex(fid, 2), 1),
          V(mesh.faceVertex(fid, 2), 2);

      double intersect_time = 0;
      bool is_intersect = IntersectTriangle(v0, v1, v2, cur_point, cur_vec,
                                            intersect_time, false);

      if (is_intersect && intersect_time >= 0) {
        next_face_id = fid;
        intersect_point = cur_point + cur_vec * intersect_time;
        next_face_orientation = mesh.tetFaceOrientation(cur_tet_id, i);
        next_face_orientation = (next_face_orientation + 1) % 2;
        next_tet_id = mesh.faceTet(next_face_id, next_face_orientation);
      }
    }
  }

  auto debug_output = [&]() {
    std::cout << "\ncur_tet_id: " << cur_tet_id << std::endl;
    std::cout << "cur_point: \n" << cur_point << std::endl;
    std::cout << "prev_vec: \n" << prev_vec << std::endl;
    std::cout << "cur_vec: \n" << cur_vec << std::endl;
    std::cout << "cur_tet_frame: \n" << cur_tet_frame << std::endl;

    bool is_pos = true;
    Eigen::Matrix3d A;
    A.row(0) = V.row(mesh.tetVertex(cur_tet_id, 1)) -
               V.row(mesh.tetVertex(cur_tet_id, 0));
    A.row(1) = V.row(mesh.tetVertex(cur_tet_id, 2)) -
               V.row(mesh.tetVertex(cur_tet_id, 0));
    A.row(2) = V.row(mesh.tetVertex(cur_tet_id, 3)) -
               V.row(mesh.tetVertex(cur_tet_id, 0));
    is_pos = A.determinant() > 0;
    std::cout << "positive orientation: " << is_pos << std::endl;

    if (s.tet_face_ids_.back() != -1) {
      int fid = mesh.tetFace(cur_tet_id, s.tet_face_ids_.back());

      int tet_id0 = mesh.faceTet(fid, 0);
      int tet_id1 = mesh.faceTet(fid, 1);
      bool is_on_bnd = (tet_id0 == -1 || tet_id1 == -1);

      Eigen::Vector3d v0, v1, v2;
      v0 << V(mesh.faceVertex(fid, 0), 0), V(mesh.faceVertex(fid, 0), 1),
          V(mesh.faceVertex(fid, 0), 2);
      v1 << V(mesh.faceVertex(fid, 1), 0), V(mesh.faceVertex(fid, 1), 1),
          V(mesh.faceVertex(fid, 1), 2);
      v2 << V(mesh.faceVertex(fid, 2), 0), V(mesh.faceVertex(fid, 2), 1),
          V(mesh.faceVertex(fid, 2), 2);

      Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
      n.normalize();

      std::cout << "boundary tet: " << is_on_bnd
                << ", dot product with pre vec: " << n.dot(prev_vec)
                << std::endl;
      std::cout << "current vec: " << n.dot(cur_vec) << std::endl;

      std::cout << "(p - v0).dot(n): " << (cur_point - v0).dot(n) << std::endl;

      Eigen::Matrix<double, 3, 2> B;
      B.col(0) = v1 - v0;
      B.col(1) = v2 - v0;

      Eigen::Matrix2d BTB = B.transpose() * B;
      Eigen::Vector2d rhs = B.transpose() * (cur_point - v0);

      Eigen::Vector2d sol = BTB.inverse() * rhs;

      std::cout << "barycentric: " << 1 - sol[0] - sol[1] << " " << sol[0]
                << " " << sol[1] << std::endl;
    }

    for (int i = 0; i < 4; i++) {
      if (i != cur_local_fid) {
        int fid = mesh.tetFace(cur_tet_id, i);
        Eigen::Vector3d v0, v1, v2;
        v0 << V(mesh.faceVertex(fid, 0), 0), V(mesh.faceVertex(fid, 0), 1),
            V(mesh.faceVertex(fid, 0), 2);
        v1 << V(mesh.faceVertex(fid, 1), 0), V(mesh.faceVertex(fid, 1), 1),
            V(mesh.faceVertex(fid, 1), 2);
        v2 << V(mesh.faceVertex(fid, 2), 0), V(mesh.faceVertex(fid, 2), 1),
            V(mesh.faceVertex(fid, 2), 2);

        double intersect_time = 0;
        bool is_intersect = IntersectTriangle(v0, v1, v2, cur_point, cur_vec,
                                              intersect_time, true);
        std::cout << "v0\n: " << v0 << "\nv1\n: " << v1 << "\nv2\n: " << v2
                  << std::endl;
        std::cout << "intersect_time: " << intersect_time
                  << ", is_intersect: " << is_intersect << std::endl;

        if (is_intersect) {
          next_face_id = fid;
          intersect_point = cur_point + cur_vec * intersect_time;
          next_face_orientation = mesh.tetFaceOrientation(cur_tet_id, i);
          next_face_orientation = (next_face_orientation + 1) % 2;
          next_tet_id = mesh.faceTet(next_face_id, next_face_orientation);
        }
      }
    }
  };

  // only for debug
  if (is_debug) {
    debug_output();
  }

  if (next_face_id == -1) {
    std::cout << "Check the code for numerical issues!" << std::endl;
    debug_output();
    return false;
  }

  if (next_tet_id == -1) {
    // reach the boundary
    return false;
  } else {
    for (int i = 0; i < 4; i++) {
      int fid = mesh.tetFace(next_tet_id, i);
      if (fid == next_face_id) {
        next_tetface_id = i;
      }
    }
  }

  s.stream_pts_.push_back(StreamPt(intersect_point, cur_vec, stream_pt_eps));
  s.tet_ids_.push_back(next_tet_id);
  s.tet_face_ids_.push_back(next_tetface_id);
  return true;
}

// Tracing streamlines
void TraceStreamlines(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const Eigen::MatrixXd& frames, int max_iter_per_trace,
                      std::vector<Streamline>& traces, int num_seeds, bool is_random_inside, double eps) {
  traces.clear();
  std::vector<std::pair<int, Eigen::Vector3d>> init_tet_ids = RandomSampleStreamPts(V, mesh, num_seeds, is_random_inside);
  int counter = 0;
  int nstreamlines = init_tet_ids.size();

  int nvecs = frames.rows() / mesh.nTets();
  std::unordered_set<StreamPt, StreamPt::HashFunction> visited_stream_pts;

  for(int vec_id = 0; vec_id < nvecs; vec_id++) {
    for (int orientation = 0; orientation < 2; orientation++) {
      int flip_sign = 1;
      if (orientation == 1) {
        flip_sign = -1;
      }
      for(int i = 0; i < nstreamlines; i++) {
        int cur_start_tet_id = init_tet_ids[i].first;

        Eigen::Vector3d cur_direction_vec = frames.row(nvecs * cur_start_tet_id + vec_id);
        cur_direction_vec *= flip_sign;

        std::cout << "processing streamline: " << counter << std::endl;
        counter++;

        Streamline t;
        t.tet_ids_.push_back(cur_start_tet_id);
        t.tet_face_ids_.push_back(-1);

        StreamPt cur_stream_pt(init_tet_ids[i].second, cur_direction_vec, eps);

        t.stream_pts_.push_back(cur_stream_pt);

        if(visited_stream_pts.count(cur_stream_pt)) {
          continue;
        }

        visited_stream_pts.insert(cur_stream_pt);

        bool is_debug = false;

        for (int j = 0; j < max_iter_per_trace; j++) {
          if (AdvanceStreamline(V, mesh, frames, t, eps, is_debug)) {
            StreamPt stream_pt = t.stream_pts_.back();
            if(visited_stream_pts.count(stream_pt)) {
              // has been visited
              t.tet_ids_.pop_back();
              t.tet_face_ids_.pop_back();
              t.stream_pts_.pop_back();
              break;
            } else {
              visited_stream_pts.insert(stream_pt);
            }
          }
          else {
            break;
          }
        }
        if (t.stream_pts_.size() >= 2) {
          // actually it is a segment
          traces.push_back(t);
        }
      }
    }
  }
}

}