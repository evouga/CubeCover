#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "SurfaceExtraction.h"
#include <Eigen/Dense>

#include "ExtractIsoSurfaces.h"
#include "polyscope/point_cloud.h"

void RenderIsoSurfaces(const std::vector<Eigen::MatrixXd>& isoVs, const std::vector<Eigen::MatrixXi>& isoFs, const std::vector<int>& iso_vals, int iso_surface_idx) {
  for(int i = 0; i < isoVs.size(); i++) {
    if(iso_surface_idx < 0 || iso_surface_idx > isoVs.size() - 1) {
      auto ps = polyscope::registerSurfaceMesh("iso surface " + std::to_string(iso_vals[i]), isoVs[i], isoFs[i]);
    } else {
      if(i == iso_surface_idx) {
        auto ps = polyscope::registerSurfaceMesh("iso surface " + std::to_string(iso_vals[i]), isoVs[i], isoFs[i]);
      } else {
        if(polyscope::hasSurfaceMesh("iso surface " + std::to_string(iso_vals[i]))) {
          polyscope::removeSurfaceMesh("iso surface " + std::to_string(iso_vals[i]));
        }
      }
    }
  }
}

std::vector<Eigen::MatrixXd> isoVs;
std::vector<Eigen::MatrixXi> isoFs;
std::vector<int> iso_vals;


Eigen::MatrixXd V;
Eigen::MatrixXi T;
CubeCover::TetMeshConnectivity mesh;
Eigen::MatrixXd values;

int iso_surface_idx = -1;


void callback() {
  ImGui::PushItemWidth(100);

  if (ImGui::CollapsingHeader("Visualization Options", ImGuiTreeNodeFlags_DefaultOpen)) {
    if(ImGui::Button("Iso Lines extraction")) {
      Eigen::MatrixXd P;
      Eigen::MatrixXi E;

      Eigen::MatrixXd P2;
      Eigen::MatrixXi E2;

      extractIsolines(V, mesh, values, P, E, P2, E2);

      auto *psCurves = polyscope::registerCurveNetwork("Isolines", P, E);
      psCurves->setRadius(0.003);
      auto *psCurves2 = polyscope::registerCurveNetwork("Bad Isolines", P2, E2);
      psCurves2->setRadius(0.003);
    }

    if(ImGui::Button("Iso Surfaces Extraction")) {
      Eigen::MatrixXd isoV;
      Eigen::MatrixXi isoF;

      CubeCover::isosurfaceSoup(V, mesh, values, isoV, isoF);

      auto *ps = polyscope::registerSurfaceMesh("Isosurfaces", isoV, isoF);
    }

    if(ImGui::Button("Multiple Surface Extraction")) {

      Eigen::MatrixXd isoV;
      Eigen::MatrixXi isoF;

      int min_val = int(values.minCoeff());
      int max_val = int(values.maxCoeff());

      for(int iso_val = min_val; iso_val < max_val; iso_val += 1) {
        CubeCover::isosurfaceSoupForSingleIsoVal(V, mesh, values, iso_val, isoV, isoF);
        isoVs.push_back(isoV);
        isoFs.push_back(isoF);
        iso_vals.push_back(iso_val);
      }

      RenderIsoSurfaces(isoVs, isoFs, iso_vals, iso_surface_idx);
    }
  }

  if (ImGui::CollapsingHeader("Iso Surface Selection Option", ImGuiTreeNodeFlags_DefaultOpen)) {
    if(ImGui::InputInt("Iso surface idx", &iso_surface_idx)) {
      RenderIsoSurfaces(isoVs, isoFs, iso_vals, iso_surface_idx);
    }
  }

  ImGui::PopItemWidth();
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> GetTetSoup(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T) {
  int ntet = T.rows();
  Eigen::MatrixXd soup_V(4 * ntet, 3);
  Eigen::MatrixXi soup_T(ntet, 4);
  for (int i = 0; i < ntet; i++) {
    for (int j = 0; j < 4; j++) {
      soup_V.row(4 * i + j) = V.row(T(i, j));
      soup_T(i, j) = 4 * i + j;
    }
  }
  return {soup_V, soup_T};
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Usage: isosurfaceviewer (.hexex file)" << std::endl;
    return -1;
  }

  std::string hexexfile = argv[1];

  if (!CubeCover::readHexEx(hexexfile, V, T, values))
  {
    std::cerr << "error reading the .hexex file" << std::endl;
    return -1;
  }

  mesh = CubeCover::TetMeshConnectivity(T);

  // make a mesh out of all of the boundary faces
  int nbdry = 0;
  int nfaces = mesh.nFaces();
  for (int i = 0; i < nfaces; i++)
  {
    if (mesh.isBoundaryFace(i))
      nbdry++;
  }
  Eigen::MatrixXi bdryF(nbdry, 3);
  int curidx = 0;
  for (int i = 0; i < nfaces; i++)
  {
    if (mesh.isBoundaryFace(i))
    {
      for (int j = 0; j < 3; j++)
      {
        bdryF(curidx, j) = mesh.faceVertex(i, j);

      }
      // fix triangle orientations
      int tet = mesh.faceTet(i, 0);
      if (tet == -1)
      {
        std::swap(bdryF(curidx, 0), bdryF(curidx, 1));
      }
      curidx++;
    }
  }

  polyscope::init();

  auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
  psMesh->setTransparency(0.2);

  Eigen::MatrixXd soup_V;
  Eigen::MatrixXi soup_T;

  std::tie(soup_V, soup_T) = GetTetSoup(V, T);

  auto tet_mesh = polyscope::registerTetMesh("tet soup mesh", soup_V, soup_T);
  std::vector<int> face_ids = {};
  for(int i = 0; i < 3; i++) {
    tet_mesh->addVertexScalarQuantity("color " + std::to_string(i),
                                      values.col(i));
    for (auto &fid : face_ids) {
      for (int j = 0; j < 4; j++) {
        if (i == 0) {
          std::cout << "color " << i << ", V" << soup_T(fid, j) << ", pos" << soup_V.row(soup_T(fid, j))
                    << ", val: " << values(soup_T(fid, j), i) << std::endl;
        }
      }
    }
  }



  // Add the callback
  polyscope::state::userCallback = callback;

  // visualize!
  polyscope::show();
}
