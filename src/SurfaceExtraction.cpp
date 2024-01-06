#include "SurfaceExtraction.h"
#include "TetMeshConnectivity.h"
#include <iostream>
#include <Eigen/Dense>

namespace CubeCover {
bool IsCrossing(double val1, double val2, double iso_val, bool flip) {
  if(!flip) {
    return val1 < iso_val && iso_val <= val2;
  } else {
    return val1 <= iso_val && iso_val < val2;
  }
}

void isosurfaceSoup(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, Eigen::MatrixXd& phivals, Eigen::MatrixXd& isoV, Eigen::MatrixXi& isoF)
{
  if (phivals.cols() != 3)
  {
    isoV.resize(0, 3);
    isoF.resize(0, 3);
    return;
  }

  int ntets = mesh.nTets();
  assert(phivals.rows() == 4 * ntets);

  std::vector<Eigen::Matrix3d> triangles;

  for (int tet = 0; tet < ntets; tet++)
  {
    for(int flip = 0; flip < 2; flip++) {
      for (int phi = 0; phi < 3; phi++)
      {
        double minphi = std::numeric_limits<double>::infinity();
        double maxphi = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < 4; i++)
        {
          double phival = phivals(4 * tet + i, phi);
          minphi = std::min(phival, minphi);
          maxphi = std::max(phival, maxphi);
        }

        int minint = int(minphi);
        int maxint = int(maxphi);
        for (int isoval = minint; isoval <= maxint; isoval++)
        {
          struct Crossing
          {
            int idx1;
            int idx2;
            Eigen::Vector3d pt;
          };
          std::vector<Crossing> crossings;
          Eigen::RowVector3d pos_pt;
          for (int idx1 = 0; idx1 < 4; idx1++)
          {
            for (int idx2 = idx1 + 1; idx2 < 4; idx2++)
            {
              int vidx1 = mesh.tetVertex(tet, idx1);
              int vidx2 = mesh.tetVertex(tet, idx2);
              double val1 = phivals(4 * tet + idx1, phi);
              double val2 = phivals(4 * tet + idx2, phi);
              if (val1 > val2)
              {
                std::swap(val1, val2);
                std::swap(vidx1, vidx2);
              }
              if (IsCrossing(val1, val2, double(isoval), flip))
              {
                double alpha = (val2 - isoval) / (val2 - val1);
                double beta = (isoval - val1) / (val2 - val1);
                Eigen::Vector3d pos = alpha * V.row(vidx1) + beta * V.row(vidx2);
                crossings.push_back({ idx1,idx2,pos });
                pos_pt = V.row(vidx2);
              }
            }
          }
          int ncrossings = crossings.size();
          assert(ncrossings == 0 || ncrossings == 3 || ncrossings == 4);

          if (ncrossings == 3)
          {
            Eigen::Matrix3d tri;
            Eigen::RowVector3d centroid(0, 0, 0);
            for (int i = 0; i < 3; i++) {
              tri.row(i) << crossings[i].pt[0], crossings[i].pt[1], crossings[i].pt[2];
              centroid += crossings[i].pt / 3;
            }
            Eigen::RowVector3d normal = (tri.row(1) - tri.row(0)).cross(tri.row(2) - tri.row(0));
            if (normal.dot(pos_pt - centroid) < 0) {
              Eigen::RowVector3d r0 = tri.row(1);
              tri.row(1) = tri.row(0);
              tri.row(0) = r0;
            }
            triangles.push_back(tri);
          }
          else if (ncrossings == 4)
          {
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
//            std::cout << centroid << std::endl;
//            std::cout << "C^T: " << centroid.transpose() << std::endl;

            for (int i = 0; i < 4; i++)
            {
              centroid += crossings[i].pt;
            }
            centroid /= 4.0;

            // fix quad orientation
            if (
                (crossings[0].idx1 != crossings[1].idx1 && crossings[0].idx1 != crossings[1].idx2)
                && (crossings[0].idx2 != crossings[1].idx1 && crossings[0].idx2 != crossings[1].idx2)
            )
            {
              std::swap(crossings[1], crossings[2]);
            }

            if (
                (crossings[1].idx1 != crossings[2].idx1 && crossings[1].idx1 != crossings[2].idx2)
                && (crossings[1].idx2 != crossings[2].idx1 && crossings[1].idx2 != crossings[2].idx2)
            )
            {
              std::swap(crossings[2], crossings[3]);
            }

            for (int offset = 0; offset < 4; offset++)
            {
              Eigen::Matrix3d tri;
              tri.row(0) << centroid[0], centroid[1], centroid[2];
              tri.row(1) << crossings[offset].pt[0], crossings[offset].pt[1], crossings[offset].pt[2];
              int op1 = (offset + 1) % 4;
              tri.row(2) << crossings[op1].pt[0], crossings[op1].pt[1], crossings[op1].pt[2];

              Eigen::RowVector3d face_center = (tri.row(0) + tri.row(1) + tri.row(2)) / 3;
              Eigen::RowVector3d face_normal = (tri.row(1) - tri.row(0)).cross(tri.row(2) - tri.row(0));

              if (face_normal.dot(pos_pt - face_center) < 0) {
                Eigen::RowVector3d r0 = tri.row(1);
                tri.row(1) = tri.row(0);
                tri.row(0) = r0;
              }

              triangles.push_back(tri);
            }
          }
        }
      }
    }

  }

  int ntris = triangles.size();
  isoV.resize(3 * ntris, 3);
  isoF.resize(ntris, 3);
  for (int i = 0; i < ntris; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      isoV.row(3 * i + j) = triangles[i].row(j);
      isoF(i, j) = 3 * i + j;
    }
  }
}

void isosurfaceSoupForSingleIsoVal(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, Eigen::MatrixXd& phivals, double isoval, Eigen::MatrixXd& isoV, Eigen::MatrixXi& isoF) {
  if (phivals.cols() != 3) {
    isoV.resize(0, 3);
    isoF.resize(0, 3);
    return;
  }

  int ntets = mesh.nTets();
  assert(phivals.rows() == 4 * ntets);

  std::vector<Eigen::Matrix3d> triangles;

  for (int tet = 0; tet < ntets; tet++) {
    for(int flip = 0; flip < 2; flip++) {
      for (int phi = 0; phi < 3; phi++) {
        struct Crossing {
          int idx1;
          int idx2;
          Eigen::Vector3d pt;
        };
        std::vector<Crossing> crossings;
        Eigen::RowVector3d pos_pt;
        for (int idx1 = 0; idx1 < 4; idx1++) {
          for (int idx2 = idx1 + 1; idx2 < 4; idx2++) {
            int vidx1 = mesh.tetVertex(tet, idx1);
            int vidx2 = mesh.tetVertex(tet, idx2);
            double val1 = phivals(4 * tet + idx1, phi);
            double val2 = phivals(4 * tet + idx2, phi);
            if (val1 > val2) {
              std::swap(val1, val2);
              std::swap(vidx1, vidx2);
            }
            if (IsCrossing(val1, val2, double(isoval), flip)) {
              double alpha = (val2 - isoval) / (val2 - val1);
              double beta = (isoval - val1) / (val2 - val1);
              Eigen::Vector3d pos =
                  alpha * V.row(vidx1) + beta * V.row(vidx2);
              crossings.push_back({idx1, idx2, pos});
              pos_pt = V.row(vidx2);
            }
          }
        }
        int ncrossings = crossings.size();
        assert(ncrossings == 0 || ncrossings == 3 || ncrossings == 4);

        if (ncrossings == 3) {
          Eigen::Matrix3d tri;
          Eigen::RowVector3d centroid(0, 0, 0);
          for (int i = 0; i < 3; i++) {
            tri.row(i) = crossings[i].pt.transpose();
            centroid += crossings[i].pt / 3;
          }
          Eigen::RowVector3d normal = (tri.row(1) - tri.row(0)).cross(tri.row(2) - tri.row(0));
          if (normal.dot(pos_pt - centroid) < 0) {
            Eigen::RowVector3d r0 = tri.row(1);
            tri.row(1) = tri.row(0);
            tri.row(0) = r0;
          }
          triangles.push_back(tri);
        } else if (ncrossings == 4) {
          Eigen::Vector3d centroid(0, 0, 0);
          for (int i = 0; i < 4; i++) {
            centroid += crossings[i].pt;
          }
          centroid /= 4.0;

          // fix quad orientation
          if ((crossings[0].idx1 != crossings[1].idx1 &&
               crossings[0].idx1 != crossings[1].idx2) &&
              (crossings[0].idx2 != crossings[1].idx1 &&
               crossings[0].idx2 != crossings[1].idx2)) {
            std::swap(crossings[1], crossings[2]);
          }

          if ((crossings[1].idx1 != crossings[2].idx1 &&
               crossings[1].idx1 != crossings[2].idx2) &&
              (crossings[1].idx2 != crossings[2].idx1 &&
               crossings[1].idx2 != crossings[2].idx2)) {
            std::swap(crossings[2], crossings[3]);
          }

          for (int offset = 0; offset < 4; offset++) {
            Eigen::Matrix3d tri;
            tri.row(0) = centroid.transpose();

            tri.row(1) = crossings[offset].pt.transpose();
            int op1 = (offset + 1) % 4;
            tri.row(2) = crossings[op1].pt.transpose();

            Eigen::RowVector3d face_center = (tri.row(0) + tri.row(1) + tri.row(2)) / 3;
            Eigen::RowVector3d face_normal = (tri.row(1) - tri.row(0)).cross(tri.row(2) - tri.row(0));

            if (face_normal.dot(pos_pt - face_center) < 0) {
              Eigen::RowVector3d r0 = tri.row(1);
              tri.row(1) = tri.row(0);
              tri.row(0) = r0;
            }

            triangles.push_back(tri);
          }
        }
      }
    }

  }

  int ntris = triangles.size();
  isoV.resize(3 * ntris, 3);
  isoF.resize(ntris, 3);
  for (int i = 0; i < ntris; i++) {
    for (int j = 0; j < 3; j++) {
      isoV.row(3 * i + j) = triangles[i].row(j);
      isoF(i, j) = 3 * i + j;
    }
  }
}
};