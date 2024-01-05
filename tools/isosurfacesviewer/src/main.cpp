#include <Eigen/Core>
#include <Eigen/Dense>
#include <filesystem>
#include <iostream>

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"

#include "TetMeshConnectivity.h"
#include "ReadFrameField.h"
#include "readMeshFixed.h"
#include "CubeCover.h"
#include "SurfaceExtraction.h"
#include "ExtractIsolines.h"
#include "TraceStreamlines.h"


enum ParametrizationType {
  kSeamless = 0,
  kIntegerGrid = 1,
};

static void RenderIsoSurfaces(const std::vector<Eigen::MatrixXd>& isoVs, const std::vector<Eigen::MatrixXi>& isoFs, const std::vector<int>& iso_vals, int iso_surface_idx) {
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

static bool IntegrateFrames(const Eigen::MatrixXd& frames, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const ParametrizationType& param_type,
                     Eigen::MatrixXi& assignments, Eigen::MatrixXd& values, double global_rescaling = 1.0) {
  CubeCover::CubeCoverOptions opt;
  if(param_type == kSeamless) {
    opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_SEAMLESS;
    opt.assignmentHandling = CubeCover::CubeCoverOptions::AssignmentHandling::AH_RECOMPUTE;
    opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FREE;

    opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_COMISO;

    opt.verbose = true;
  } else {
    opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_INTEGERGRID;
    opt.assignmentHandling = CubeCover::CubeCoverOptions::AssignmentHandling::AH_RECOMPUTE;
    opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER;
    opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_COMISO;

    // set to something non-zero if you want curl-correction. 1.0 == 100% change in the input frames allowed.
    opt.curlCorrection = 0.0;

    opt.verbose = true;
  }
  Eigen::MatrixXd to_round = global_rescaling * frames;
  return CubeCover::cubeCover(V, T, to_round, assignments, values, opt);
}

static void hsv_to_rgb(
    const double & h, const double& s, const double& v,
    double& r, double& g, double & b)
{
  // From medit
  double f,p,q,t,hh;
  int    i;
  // shift the hue to the range [0, 360] before performing calculations
  hh = ((360 + ((int)h % 360)) % 360) / 60.;
  i = (int)std::floor(hh);    /* largest int <= h     */
  f = hh - i;                    /* fractional part of h */
  p = v * (1.0 - s);
  q = v * (1.0 - (s * f));
  t = v * (1.0 - (s * (1.0 - f)));

  switch(i) {
  case 0: r = v; g = t; b = p; break;
  case 1: r = q; g = v; b = p; break;
  case 2: r = p; g = v; b = t; break;
  case 3: r = p; g = q; b = v; break;
  case 4: r = t; g = p; b = v; break;
  case 5: r = v; g = p; b = q; break;
  }
}

static Eigen::MatrixXd PaintPhi(const Eigen::VectorXd& phi, Eigen::VectorXd* brightness = nullptr)      // brightness should between 0 and 1
{
  int nverts = phi.size();
  Eigen::MatrixXd color(nverts, 3);
  for (int i = 0; i < nverts; i++)
  {
    double r, g, b;
    //            double h = 360.0 * phi[i] / 2.0 / M_PI + 120;
    double h = 360.0 * phi[i] / 2.0 / M_PI;
    h = 360 + ((int)h % 360); // fix for libigl bug
    double s = 1.0;
    double v = 0.5;
    if(brightness)
    {
      double r = (*brightness)(i);
      v = r * r / (r * r + 1);
    }
    //                v = (*brightness)(i);
    hsv_to_rgb(h, s, v, r, g, b);
    color(i, 0) = r;
    color(i, 1) = g;
    color(i, 2) = b;
  }
  return color;
}

static void RenderScalarFields(polyscope::VolumeMesh* tet_mesh, const Eigen::MatrixXd& values) {
  for(int i = 0; i < 3; i++) {
    Eigen::MatrixXd vertex_color = PaintPhi(values.col(i));
    tet_mesh->addVertexColorQuantity("color " + std::to_string(i),
                                      vertex_color);
  }
}

static void RenderStreamlines(const std::vector<Streamline>& traces, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T) {
  Eigen::VectorXd tet_colors;
  tet_colors.resize(T.rows());
  tet_colors.setZero();
  int ntets_mesh = T.rows();
  int ntraces = traces.size();

  std::vector<Eigen::Vector4i> streamline_tets;


  for (int tid = 0; tid < ntraces; tid++)
  {
    int nsteps = traces.at(tid).tetIds.size();
    for (int i = 0; i < nsteps; i++ )
    {
      auto cur_trace = traces.at(tid);
      int cur_tet_id = cur_trace.tetIds.at(i);
      streamline_tets.push_back(T.row(cur_tet_id));
      // std::cout << "cur_tet_id: " << cur_tet_id << std::endl;
      tet_colors(cur_tet_id) = 1.;
    }
  }

  Eigen::MatrixXd cur_points;
  Eigen::MatrixXd points;
  int nfam = 6;

  for (int cur_fam_id = 0; cur_fam_id < nfam; cur_fam_id++)
  {
    int iter = 0;
    std::vector<Eigen::Vector2d> cur_line_edges;
    std::vector<Eigen::Vector3d> cur_points;

    for (int tid = 0; tid < ntraces; tid++)
    {
      int tid_fam_id = tid % nfam;
      if (tid_fam_id == cur_fam_id)
      {
        //       streamline_tets.clear();
        int nsteps = traces.at(tid).points.size();

        Eigen::Vector3d first_edge = traces.at(tid).points.at(1) - traces.at(tid).points.at(0);
        double fen = first_edge.norm();

        int cur_len = 0;
        bool addLast = true;
        for (int i = 0; i < nsteps-1; i++ )
        {
          Eigen::Vector3d edge = traces.at(tid).points.at(i) - traces.at(tid).points.at(i+1);
          if (edge.norm() > fen * 5 )
          {
            addLast = false;
            break;
          }
          cur_points.push_back( traces.at(tid).points.at(i) );
          cur_len++;

        }
        if (addLast)
        {
          cur_points.push_back( traces.at(tid).points.at(nsteps-1) );
          cur_len++;
        }




        for (int i = 0; i < cur_len-1; i++ )
        {
          cur_line_edges.push_back(Eigen::Vector2d(iter+i, iter+i+1) );
        }
        iter = iter + cur_len;

      }



    }

    auto *single_streamline = polyscope::registerCurveNetwork("streamline" + std::to_string(cur_fam_id), cur_points, cur_line_edges);
    single_streamline->setTransparency(1);
    single_streamline->setRadius(0.003);


  }
}

std::vector<Eigen::MatrixXd> ExtractFrameVectors(int ntets, const Eigen::MatrixXd& frames) {
  int nvecs = frames.rows() / ntets;
  std::vector<Eigen::MatrixXd> to_ret;
  for(int i = 0; i < nvecs; i++) {
    Eigen::MatrixXd frame_i(ntets, 3);
    for(int j = 0; j < ntets; j++) {
      frame_i.row(j) = frames.row(j * nvecs + i);
    }
    to_ret.emplace_back(frame_i);
  }
  return to_ret;
}

std::vector<Eigen::MatrixXd> isoVs;
std::vector<Eigen::MatrixXi> isoFs;
std::vector<int> iso_vals;


Eigen::MatrixXd V;
Eigen::MatrixXi T;
CubeCover::TetMeshConnectivity mesh;

Eigen::MatrixXd frames, frames_to_trace;
Eigen::MatrixXi assignments;
Eigen::MatrixXd values;

double global_rescaling = 1.0;

int iso_surface_idx = -1;

ParametrizationType param_type = kSeamless;

std::string save_folder;


void callback() {
  ImGui::PushItemWidth(100);

  if (ImGui::CollapsingHeader("Integration Options", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Combo("Paramaterization Type", (int*)&param_type, "Seamless\0Integer grid\0");
    ImGui::InputDouble("Global Rescale", &global_rescaling);
    if(ImGui::Button("Integrate Frames")) {
      if(!IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling)) {
        std::cout << "cube cover failed!" << std::endl;
      } else {
        auto tet_mesh = polyscope::getVolumeMesh("tet soup mesh");
        RenderScalarFields(tet_mesh, values);
      }
    }
  }

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
    if(ImGui::Button("Trace Stream lines")) {
      std::vector<Streamline> traces;
      int max_iter_per_trace = 700;

      Eigen::VectorXi init_tet_ids;
      CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames_to_trace, assignments, true);

      field->computeLocalAssignments();
      field->combAssignments();

      traceStreamlines(V, mesh, *field, init_tet_ids, max_iter_per_trace, traces);

      RenderStreamlines(traces, V, T);
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
  if (argc != 3 && argc != 4)
  {
    std::cerr << "Usage: isosurfaceviewer (.mesh file) (.bfra file or .fra file) [save folder]" << std::endl;
    return -1;
  }

  std::string meshfile = argv[1];
  std::string frafile = argv[2];
  std::string permfile;

  std::string outfile;
  if (argc == 4)
    save_folder = argv[3];

  Eigen::MatrixXi F;
  if (!CubeCover::readMESH(meshfile, V, T, F))
  {
    std::cerr << "could not read .mesh file " << meshfile << std::endl;
    return -1;
  }

  if (!CubeCover::readFrameField(frafile, permfile, T, frames, assignments, true))
  {
    std::cerr << "could not read frames/permutations" << std::endl;
    return -1;
  }

  int vpe = frames.rows() / T.rows();
  if(vpe != 3) {
    Eigen::MatrixXd ext_frames(3 * T.rows(), 3);

    for(int i = 0; i < T.rows(); i++) {
      if(vpe == 0) {
        ext_frames.row(3 * i + 0) << 1, 0, 0;
        ext_frames.row(3 * i + 1) << 0, 1, 0;
        ext_frames.row(3 * i + 2) << 0, 0, 1;
      } else if (vpe == 1) {
        ext_frames.row(3 * i + 0) << frames.row(vpe * i + 0);
        ext_frames.row(3 * i + 1) << 0, 1, 0;
        ext_frames.row(3 * i + 2) << 0, 0, 1;
      } else if (vpe == 2) {
        ext_frames.row(3 * i + 0) << frames.row(vpe * i + 0);
        ext_frames.row(3 * i + 1) << frames.row(vpe * i + 1);
        ext_frames.row(3 * i + 2) << 0, 0, 1;
      } else {
        for(int j = 0; j < 3; j++) {
          ext_frames.row(3 * i + j) = frames.row(vpe * i + j);
        }
      }
    }

    frames = std::move(ext_frames);
  }

  frames_to_trace = frames;
  for(int i = 0; i < T.rows(); i++) {
    for(int j = 0; j < 3; j++) {
      Eigen::Vector3d v1 = frames.row(3 * i + (j + 1) % 3);
      Eigen::Vector3d v2 = frames.row(3 * i + (j + 2) % 3);
      Eigen::Vector3d v0 = v1.cross(v2);
      frames_to_trace.row(3 * i + j) = v0;
    }
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

  IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling);

  std::vector<Eigen::MatrixXd> frame_vec = ExtractFrameVectors(T.rows(), frames);

  polyscope::init();

  auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
  psMesh->setTransparency(0.2);
  psMesh->setEnabled(false);

  Eigen::MatrixXd soup_V;
  Eigen::MatrixXi soup_T;

  std::tie(soup_V, soup_T) = GetTetSoup(V, T);

  auto tet_mesh = polyscope::registerTetMesh("tet soup mesh", soup_V, soup_T);
  std::vector<int> face_ids = {};
  for(int i = 0; i < 3; i++) {
    tet_mesh->addCellVectorQuantity("frame " + std::to_string(i), frame_vec[i]);
  }

  RenderScalarFields(tet_mesh, values);

  // Add the callback
  polyscope::state::userCallback = callback;

  // visualize!
  polyscope::show();
}
