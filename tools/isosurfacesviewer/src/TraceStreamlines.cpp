#include "TraceStreamlines.h"
#include "TetMeshConnectivity.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

#include "FrameField.h"

#include <set>


static bool intersectPlane(const Eigen::Vector3d &n, const Eigen::Vector3d &p0, const Eigen::Vector3d &l0, const Eigen::Vector3d &l, double &t)
{
  // assuming vectors are all normalized
  double denom = n.dot(l);
  if (denom > 1e-6) {
    Eigen::Vector3d p0l0 = p0 - l0;
    t = p0l0.dot(n) / denom;
    //    t = p0l0.dot(n) / denom - .001; // in case tet is degenerate, don't quite hit the surface
    return (t >= 0);
  }

  return false;
}

// check whether a line: l(t) = p0 + t * dir is intersection with the triangle v0, v1, v2;
// If intersection happens at p = p0 + t * dir, we have p - v0 = c0 (v1 - v0) + c1 (v2- v0) => c0 * (v1 - v0)  + c1 * (v2 - v0) - t * dir = p0 - v0
// We can solve the linear system [v1 - v0, v2 - v0, -dir] [c0, c1, t]^T = [p0 - v0]. This can be solved if dir is not on the triangle plane
static bool intersectTriangle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& p0, const Eigen::Vector3d& dir, double& t, bool is_debug = false) {
  Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
  n.normalize();

  // parallel case
  if(std::abs(n.dot(dir)) < 1e-6) {
    return false;
  }

  Eigen::Matrix3d A;
  A.col(0) = v1 - v0;
  A.col(1) = v2 - v0;
  A.col(2) = -dir;

  Eigen::Vector3d rhs = p0 - v0;

  Eigen::Vector3d sol = A.inverse() * rhs;

  t = sol[2];
  return (sol[0] >= 0) && (sol[1] >= 0) && (sol[0] + sol[1]) <= 1;
}


static bool advanceStreamline( const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const CubeCover::FrameField& frames,
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

  // only for debug
  if(is_debug) {
    std::cout << "\ncur_tet_id: " << cur_tet_id << std::endl;
    std::cout << "cur_point: " << cur_point.transpose() << std::endl;
    std::cout << "prev_vec: " << prev_vec.transpose() << std::endl;
    std::cout << "cur_vec: " << cur_vec.transpose() << std::endl;
    std::cout << "cur_tet_frame: \n" << cur_tet_frame << std::endl;

    int fid = mesh.tetFace(cur_tet_id, s.tetFaceIds.back());

    int tet_id0 = mesh.faceTet(fid, 0);
    int tet_id1 = mesh.faceTet(fid, 1);
    bool is_on_bnd = (tet_id0 == -1 || tet_id1 == -1);

    Eigen::Vector3d v0 = V.row(mesh.faceVertex(fid, 0));
    Eigen::Vector3d v1 = V.row(mesh.faceVertex(fid, 1));
    Eigen::Vector3d v2 = V.row(mesh.faceVertex(fid, 2));

    Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
    n.normalize();

    std::cout << "boundary tet: " << is_on_bnd << ", dot product with pre vec: " << n.dot(prev_vec) << std::endl;
    std::cout << "current vec: " << n.dot(cur_vec) << std::endl;

    bool is_pos = true;
    Eigen::Matrix3d A;
    A.row(0) = V.row(mesh.tetVertex(cur_tet_id, 1)) - V.row(mesh.tetVertex(cur_tet_id, 0));
    A.row(1) = V.row(mesh.tetVertex(cur_tet_id, 2)) - V.row(mesh.tetVertex(cur_tet_id, 0));
    A.row(2) = V.row(mesh.tetVertex(cur_tet_id, 3)) - V.row(mesh.tetVertex(cur_tet_id, 0));
    is_pos = A.determinant() > 0;
    std::cout << "positive orientation: " << is_pos << std::endl;
    std::cout << "(p - v0).dot(n): " << (cur_point - v0).dot(n) << std::endl;

    Eigen::Matrix<double, 3, 2> B;
    B.col(0) = v1 - v0;
    B.col(1) = v2 - v0;

    Eigen::Matrix2d BTB = B.transpose() * B;
    Eigen::Vector2d rhs = B.transpose() * (cur_point - v0);

    Eigen::Vector2d sol = BTB.inverse() * rhs;

    std::cout << "barycentric: " << 1 - sol[0] - sol[1] << " " << sol.transpose() << std::endl;
  }

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
      Eigen::Vector3d v0 = V.row(mesh.faceVertex(fid, 0));
      Eigen::Vector3d v1 = V.row(mesh.faceVertex(fid, 1));
      Eigen::Vector3d v2 = V.row(mesh.faceVertex(fid, 2));

      double intersect_time = 0;
      bool is_intersect = intersectTriangle(v0, v1, v2, cur_point, cur_vec, intersect_time, is_debug);

      if(is_debug) {
        std::cout << "intersect_time: " << intersect_time << ", is_intersect: " << is_intersect << std::endl;
      }

      if(is_intersect) {
        next_face_id = fid;
        intersect_point = cur_point + cur_vec*intersect_time;
        next_face_orientation = mesh.tetFaceOrientation(cur_tet_id, i);
        next_face_orientation = (next_face_orientation + 1) % 2;
        next_tet_id = mesh.faceTet(next_face_id, next_face_orientation);
      }
    }
  }

  if ( next_face_id == -1)
  {
    std::cout << "Something wierd happened! Ray already at the boundary!" << std::endl;
    std::cout << "current tet id: " << cur_tet_id << std::endl;

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

static Eigen::VectorXi randomSample(int ntets, int nsamples) {
  if(ntets <= nsamples) {
    Eigen::VectorXi sample_ids(ntets);
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

    Eigen::VectorXi to_ret(tmp_vec.size());
    for(int i=0; i < tmp_vec.size(); i++) {
      to_ret[i] = tmp_vec[i];
    }
    return to_ret;
  }
}

static Eigen::VectorXi sequentialSample(const CubeCover::TetMeshConnectivity& mesh, int ntets, int nsamples) {
  std::vector<int> sampledTets;
  std::unordered_set<int> visitedTets;

  for (int i = 0; i < std::min(nsamples, ntets); i++) {
    auto it = visitedTets.find(i);
    if (it == visitedTets.end()) {
      visitedTets.insert(i);
      sampledTets.push_back(i);
      // insert the neighbors.
      for(int j = 0; j < 4; j++) {
        int f1 = mesh.faceTet(mesh.tetFace(i, j), 0);
        int f2 = mesh.faceTet(mesh.tetFace(i, j), 1);

        if(f1 != i && f1 != -1) {
          visitedTets.insert(f1);
        } else if (f2 != i && f2 != -1) {
          visitedTets.insert(f2);
        }
      }
    }
  }
  int samplesize = sampledTets.size();
  Eigen::VectorXi to_ret(samplesize);
  for(int k = 0; k < samplesize; k++) {
    to_ret(k) = sampledTets.at(k);
  }
  return to_ret;
}


void traceStreamlines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const CubeCover::FrameField& frames,
                      Eigen::VectorXi& init_tet_ids,
                      int max_iter_per_trace,
                      std::vector<Streamline>& traces)
{
  traces.clear();

  if (init_tet_ids.size() < 1)
  {
    std::vector<int> sampledTets;
    std::unordered_set<int> visitedTets;

    int ntets = mesh.nTets();
//    init_tet_ids = sequentialSample(mesh, ntets, 100);
    init_tet_ids = randomSample(ntets, 100);
  }

  int counter = 0;

  int nstreamlines = init_tet_ids.rows();
  for (int i = 0; i < nstreamlines; i++)
  {
    int cur_start_tet_id = init_tet_ids(i);
    Eigen::MatrixXd cur_tet_frame = frames.tetFrame(cur_start_tet_id);
    int nvecs = cur_tet_frame.rows();

    for (int vec_id = 0; vec_id < nvecs; vec_id++)
    {
      Eigen::Vector3d cur_direction_vec = cur_tet_frame.row(vec_id);

      for (int orientation = 0; orientation < 2; orientation++)
      {
        if (orientation == 1)
        {
          cur_direction_vec = -cur_direction_vec;
        }

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

        // Eigen::Vector3d cur_vec = cur_tet_frame.row(0);

        bool canContinue = true;

        for (int j = 0; j < max_iter_per_trace; j++)
        {
          if (canContinue)
          {
            canContinue = advanceStreamline(V, mesh, frames, t, false);
          }
          else
          {
            break;
          }
        }
        if (t.tetIds.size() >= 1)
        {
          traces.push_back(t);
        }

      }
    }
  }
}