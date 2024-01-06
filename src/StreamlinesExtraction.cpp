#include "StreamlinesExtraction.h"

#include <cassert>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>



namespace CubeCover {
// check whether a line: l(t) = p0 + t * dir is intersection with the triangle v0, v1, v2;
// If intersection happens at p = p0 + t * dir, we have p - v0 = c0 (v1 - v0) + c1 (v2- v0) => c0 * (v1 - v0)  + c1 * (v2 - v0) - t * dir = p0 - v0
// We can solve the linear system [v1 - v0, v2 - v0, -dir] [c0, c1, t]^T = [p0 - v0]. This can be solved if dir is not on the triangle plane
static bool IntersectTriangle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& p0, const Eigen::Vector3d& dir, double& t, bool is_debug = false) {
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


static bool AdvanceStreamline( const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const FrameField& frames,
                              Streamline& s, bool is_debug = false)
{
  int cur_tet_id = s.tetIds.back();

  Eigen::MatrixXd cur_tet_frame = frames.tetFrame(cur_tet_id);
  Eigen::Vector3d cur_point = s.points.back();
  Eigen::Vector3d prev_vec = s.directionVecs.back();

  double max_match = 0;
  Eigen::Vector3d first_frame_vec = cur_tet_frame.row(0);
  Eigen::Vector3d cur_vec = prev_vec.cross(first_frame_vec);
  for (int i = 0; i < cur_tet_frame.rows(); i++)
  {
    Eigen::Vector3d cur_dir = cur_tet_frame.row(i);
    double cur_match = cur_dir.dot(prev_vec) / cur_dir.norm() / prev_vec.norm();

    if (std::abs(cur_match) > max_match)
    {
      max_match = std::abs(cur_match);
      cur_vec = cur_tet_frame.row(i);
      if (cur_match < 0)
        cur_vec = -cur_vec;
    }
  }
  cur_vec.normalize();

  int cur_local_fid = s.tetFaceIds.back();

  int next_tetface_id = -1;
  int next_face_id = -1;
  int next_tet_id = -1;
  Eigen::Vector3d intersect_point = cur_point;
  int next_face_orientation = -1;

  for (int i = 0; i < 4; i++)
  {
    if (i != cur_local_fid)
    {
      int fid = mesh.tetFace(cur_tet_id, i);
      Eigen::Vector3d v0, v1, v2;
      v0 << V(mesh.faceVertex(fid, 0), 0), V(mesh.faceVertex(fid, 0), 1), V(mesh.faceVertex(fid, 0), 2);
      v1 << V(mesh.faceVertex(fid, 1), 0), V(mesh.faceVertex(fid, 1), 1), V(mesh.faceVertex(fid, 1), 2);
      v2 << V(mesh.faceVertex(fid, 2), 0), V(mesh.faceVertex(fid, 2), 1), V(mesh.faceVertex(fid, 2), 2);

      double intersect_time = 0;
      bool is_intersect = IntersectTriangle(v0, v1, v2, cur_point, cur_vec, intersect_time, false);

      if(is_intersect) {
        next_face_id = fid;
        intersect_point = cur_point + cur_vec*intersect_time;
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
    A.row(0) = V.row(mesh.tetVertex(cur_tet_id, 1)) - V.row(mesh.tetVertex(cur_tet_id, 0));
    A.row(1) = V.row(mesh.tetVertex(cur_tet_id, 2)) - V.row(mesh.tetVertex(cur_tet_id, 0));
    A.row(2) = V.row(mesh.tetVertex(cur_tet_id, 3)) - V.row(mesh.tetVertex(cur_tet_id, 0));
    is_pos = A.determinant() > 0;
    std::cout << "positive orientation: " << is_pos << std::endl;

    if(s.tetFaceIds.back() != -1) {
      int fid = mesh.tetFace(cur_tet_id, s.tetFaceIds.back());

      int tet_id0 = mesh.faceTet(fid, 0);
      int tet_id1 = mesh.faceTet(fid, 1);
      bool is_on_bnd = (tet_id0 == -1 || tet_id1 == -1);

      Eigen::Vector3d v0, v1, v2;
      v0 << V(mesh.faceVertex(fid, 0), 0), V(mesh.faceVertex(fid, 0), 1), V(mesh.faceVertex(fid, 0), 2);
      v1 << V(mesh.faceVertex(fid, 1), 0), V(mesh.faceVertex(fid, 1), 1), V(mesh.faceVertex(fid, 1), 2);
      v2 << V(mesh.faceVertex(fid, 2), 0), V(mesh.faceVertex(fid, 2), 1), V(mesh.faceVertex(fid, 2), 2);

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

      std::cout << "barycentric: " << 1 - sol[0] - sol[1] << " "
                << sol[0] << " " << sol[1] << std::endl;
    }

    for (int i = 0; i < 4; i++)
    {
      if (i != cur_local_fid)
      {
        int fid = mesh.tetFace(cur_tet_id, i);
        Eigen::Vector3d v0, v1, v2;
        v0 << V(mesh.faceVertex(fid, 0), 0), V(mesh.faceVertex(fid, 0), 1), V(mesh.faceVertex(fid, 0), 2);
        v1 << V(mesh.faceVertex(fid, 1), 0), V(mesh.faceVertex(fid, 1), 1), V(mesh.faceVertex(fid, 1), 2);
        v2 << V(mesh.faceVertex(fid, 2), 0), V(mesh.faceVertex(fid, 2), 1), V(mesh.faceVertex(fid, 2), 2);

        double intersect_time = 0;
        bool is_intersect = IntersectTriangle(v0, v1, v2, cur_point, cur_vec, intersect_time, true);
        std::cout << "v0\n: " << v0 << "\nv1\n: " << v1 << "\nv2\n: " << v2 << std::endl;
        std::cout << "intersect_time: " << intersect_time << ", is_intersect: " << is_intersect << std::endl;

        if(is_intersect) {
          next_face_id = fid;
          intersect_point = cur_point + cur_vec*intersect_time;
          next_face_orientation = mesh.tetFaceOrientation(cur_tet_id, i);
          next_face_orientation = (next_face_orientation + 1) % 2;
          next_tet_id = mesh.faceTet(next_face_id, next_face_orientation);
        }
      }
    }
  };

  // only for debug
  if(is_debug) {
    debug_output();
  }

  if ( next_face_id == -1)
  {
    std::cout << "Check the code for numerical issues!" << std::endl;
    debug_output();

    s.directionVecs.pop_back();
    s.points.pop_back();
    s.faceIds.pop_back();
    return false;
  }

  s.directionVecs.push_back(cur_vec);
  s.points.push_back(intersect_point);
  s.faceIds.push_back(next_face_id);

  if ( next_tet_id == -1)
  {
    return false;
  }
  else
  {
    for (int i = 0; i < 4; i++)
    {
      int fid = mesh.tetFace(next_tet_id, i);
      if (fid == next_face_id)
      {
        next_tetface_id = i;
      }
    }
  }

  s.tetIds.push_back(next_tet_id);
  s.tetFaceIds.push_back(next_tetface_id);
  return true;

}

static std::vector<int> RandomSample(int ntets, int nsamples) {
  if(ntets <= nsamples) {
    std::vector<int> sample_ids(ntets);
    for(int i = 0; i < ntets; i++) {
      sample_ids[i] = i;
    }
    return sample_ids;
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

    return tmp_vec;
  }
}

static std::vector<int> SequentialSample(const TetMeshConnectivity& mesh, int ntets, int nsamples) {
  std::vector<int> sampled_tets;
  std::unordered_set<int> visited_tets;

  for (int i = 0; i < std::min(nsamples, ntets); i++) {
    auto it = visited_tets.find(i);
    if (it == visited_tets.end()) {
      visited_tets.insert(i);
      sampled_tets.push_back(i);
      // insert the neighbors.
      for(int j = 0; j < 4; j++) {
        int f1 = mesh.faceTet(mesh.tetFace(i, j), 0);
        int f2 = mesh.faceTet(mesh.tetFace(i, j), 1);

        if(f1 != i && f1 != -1) {
          visited_tets.insert(f1);
        } else if (f2 != i && f2 != -1) {
          visited_tets.insert(f2);
        }
      }
    }
  }
  return sampled_tets;
}

static bool IsCrossing(double val1, double val2, double iso_val, bool flip) {
  if(!flip) {
    return val1 < iso_val && iso_val <= val2;
  } else {
    return val1 <= iso_val && iso_val < val2;
  }
}

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

void TraceStreamlines(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, const FrameField& frames,
                      std::vector<int>& init_tet_ids,
                      int max_iter_per_trace,
                      std::vector<Streamline>& traces) {
  traces.clear();
  if (init_tet_ids.empty()) {
    init_tet_ids = RandomSample(mesh.nTets(), 100);
  }
  int counter = 0;
  int nstreamlines = init_tet_ids.size();

  struct TetDirect {
    int id_;
    Eigen::Vector3d dir_;

    TetDirect(int id, const Eigen::Vector3d& dir) : id_(id), dir_(dir) {
      dir_.normalize();
    }

    bool operator==(const TetDirect& other) const {
      return this->id_ == other.id_ && (this->dir_ - other.dir_).norm() == 0;
    }

    struct HashFunction {
      size_t operator()(const TetDirect& tet_dir) const {
        std::size_t h1 = std::hash<int>{}(tet_dir.id_);
        std::size_t h2 = std::hash<double>{}(tet_dir.dir_[0]);
        std::size_t h3 = std::hash<double>{}(tet_dir.dir_[1]);
        std::size_t h4 = std::hash<double>{}(tet_dir.dir_[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
      }
    };
  };


  int nvecs = frames.tetFrame(0).rows();
  std::unordered_set<TetDirect, TetDirect::HashFunction> visited_tet_ids;

  for(int vec_id = 0; vec_id < nvecs; vec_id++) {
    for (int orientation = 0; orientation < 2; orientation++) {
      bool flip_sign = 1;
      if (orientation == 1) {
        flip_sign = -1;
      }
      for(int i = 0; i < nstreamlines; i++) {
        int cur_start_tet_id = init_tet_ids[i];
        Eigen::Vector3d cur_direction_vec = frames.tetFrame(cur_start_tet_id).row(vec_id);
        cur_direction_vec *= flip_sign;

        TetDirect cur_tet_dir(cur_start_tet_id, cur_direction_vec);
//        if(visited_tet_ids.count(cur_tet_dir)) {
//          continue;
//        }

        std::cout << "processing streamline: " << counter << std::endl;
        counter++;

        Streamline t;
        t.tetIds.push_back(cur_start_tet_id);
        t.tetFaceIds.push_back(-1);

        int cur_face_id = -1;
        t.faceIds.push_back(cur_face_id);
        t.directionVecs.push_back(cur_direction_vec);


        Eigen::Vector4i cur_tet;
        cur_tet(0) = mesh.tetVertex(cur_start_tet_id, 0);
        cur_tet(1) = mesh.tetVertex(cur_start_tet_id, 1);
        cur_tet(2) = mesh.tetVertex(cur_start_tet_id, 2);
        cur_tet(3) = mesh.tetVertex(cur_start_tet_id, 3);

        Eigen::Vector3d cur_point;
        Eigen::Vector3d v0 = V.row(cur_tet(0));
        Eigen::Vector3d v1 = V.row(cur_tet(1));
        Eigen::Vector3d v2 = V.row(cur_tet(2));
        Eigen::Vector3d v3 = V.row(cur_tet(3));
        cur_point = ( v0 + v1 + v2 + v3 ) / 4.;
        t.points.push_back(cur_point);

        visited_tet_ids.insert(cur_tet_dir);

        bool canContinue = true;

        bool is_debug = false;

        for (int j = 0; j < max_iter_per_trace; j++) {
          if (canContinue) {
            canContinue = AdvanceStreamline(V, mesh, frames, t, is_debug);  // during the AdvanceStreamline we check the -dir and dir for best match
            TetDirect tet_dir(t.tetIds.back(), t.directionVecs.back());
            if(visited_tet_ids.count(tet_dir)) {
              // has been visited
              canContinue = false;
              t.tetIds.pop_back();
              t.directionVecs.pop_back();
              t.points.pop_back();
              t.faceIds.pop_back();
              t.tetFaceIds.pop_back();
              break;
            } else {
              visited_tet_ids.insert(tet_dir);
            }
          }
          else {
            break;
          }
        }
        if (t.points.size() >= 2) {
          traces.push_back(t);
        }

      }
    }

  }
}
}