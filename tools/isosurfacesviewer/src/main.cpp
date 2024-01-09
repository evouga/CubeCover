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
#include "ReadHexEx.h"
#include "CubeCover.h"
#include "SurfaceExtraction.h"
#include "StreamlinesExtraction.h"
#include "IsolinesExtraction.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950
#endif // !M_PI



enum ParametrizationType {
  kSeamless = 0,
  kIntegerGrid = 1,
};

enum StreamLineType {
  kInit = 0,
  kGradient = 1,
  kInitPerp = 2,
  kInitBestMatchGrad = 3,
};

enum StreamLineTracingType {
    kRandomCentroid = 0,    // randomly sample the tet ids and use its centroid as the starting point
    kRandom = 1,            // randomly sample the tet ids and use a random point inside tets for starting point
    kGridPt = 2,            // use integer grid points
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

static Eigen::MatrixXd ComputeGradient(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values) {
  int nvecs = values.cols();
  int ntets = mesh.nTets();

  // assume that the values are defined on the tet soup
  auto get_soup_vid = [&](int tet_id, int j) {
    return 4 * tet_id + j;    // this is corresponding to the tet(tet_id, j)
  };

  Eigen::MatrixXd grad_fields(ntets * nvecs, 3);

  for(int vec_id = 0; vec_id < nvecs; vec_id++) {
    for(int tet_id = 0; tet_id < ntets; tet_id++) {
      Eigen::Matrix3d A;
      Eigen::Vector3d rhs;
      for(int j = 0; j < 3; j++) {
        A.row(j) = V.row(mesh.tetVertex(tet_id, j + 1)) - V.row(mesh.tetVertex(tet_id, 0));
        rhs[j] = values(get_soup_vid(tet_id, j + 1), vec_id) - values(get_soup_vid(tet_id, 0), vec_id);
      }
      if(std::abs(A.determinant()) > 1e-10) {
        Eigen::Vector3d grad = A.inverse() * rhs;
        grad_fields.row(nvecs * tet_id + vec_id) << grad[0], grad[1], grad[2];
      } else {
        grad_fields.row(nvecs * tet_id + vec_id).setZero();
      }
    }
  }
  return grad_fields;
}

std::vector<Eigen::VectorXd> GetFrameDifference(const Eigen::MatrixXd& frame0, const Eigen::MatrixXd& frame1, int nvecs) {
  int neles = frame0.rows() / nvecs;
  std::vector<Eigen::VectorXd> errs(nvecs, Eigen::VectorXd::Zero(neles));
  for(int vec_id = 0; vec_id < nvecs; vec_id++) {
    for(int ele_id = 0; ele_id < neles; ele_id++) {
      Eigen::RowVector3d v0 = frame0.row(nvecs * ele_id + vec_id);

      double min_err = std::numeric_limits<double>::max();
      for(int vec_id1 = 0; vec_id1 < nvecs; vec_id1++) {
        Eigen::RowVector3d v1 = frame1.row(nvecs * ele_id + vec_id1);
        double err = (v0 - v1).norm();
        min_err = std::min(min_err, err);

        err = (v0 + v1).norm();
        min_err = std::min(min_err, err);
      }

      errs[vec_id][ele_id] = min_err;
    }
  }
  return errs;
}

std::vector<Eigen::MatrixXd> GetBestMatchFrames(const std::vector<Eigen::MatrixXd>& frame_list, const Eigen::MatrixXd& frame1) {
  if(frame_list.empty()) {
    return {};
  }

  int neles = frame_list[0].rows();
  int nvecs = frame_list.size();
  std::vector<Eigen::MatrixXd> best_frame_list_1(nvecs, Eigen::MatrixXd::Zero(neles, 3));
  for(int vec_id = 0; vec_id < nvecs; vec_id++) {
    for(int ele_id = 0; ele_id < neles; ele_id++) {
      Eigen::RowVector3d v0 = frame_list[vec_id].row(ele_id);

      double min_err = std::numeric_limits<double>::max();
      for(int vec_id1 = 0; vec_id1 < nvecs; vec_id1++) {
        Eigen::RowVector3d v1 = frame1.row(nvecs * ele_id + vec_id1);
        double err = (v0 - v1).norm();
        if(err < min_err) {
          best_frame_list_1[vec_id].row(ele_id) = v1;
          min_err = err;
        }
        err = (v0 + v1).norm();
        if(err < min_err) {
          best_frame_list_1[vec_id].row(ele_id) = -v1;
          min_err = err;
        }
      }
    }
  }
  return best_frame_list_1;
}

static void RenderScalarFields(polyscope::VolumeMesh* tet_mesh, const Eigen::MatrixXd& values) {
  for(int i = 0; i < 3; i++) {
    Eigen::MatrixXd vertex_color = PaintPhi(values.col(i));
    tet_mesh->addVertexColorQuantity("color " + std::to_string(i),
                                      vertex_color);
  }
}

static Eigen::Vector3d SegmentColor(const Eigen::Vector3d& dir, Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity()) {
  Eigen::Vector3d dir_hat = dir.normalized();
  dir_hat = rot_mat * dir_hat;
  double phi = atan2(dir_hat[1], dir_hat[0]) + M_PI;
  double theta = acos(dir_hat[1]);

  double h = 360 / 2 / M_PI * phi;
  double v = 1;
  double s = std::sin(theta);

  Eigen::Vector3d rgb;
  hsv_to_rgb(h, s, v, rgb[0], rgb[1], rgb[2]);
  return rgb;
}

static void RenderStreamlines(const std::vector<CubeCover::Streamline>& traces, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, int nframes = 1, std::string name = "", Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity()) {
  if(traces.empty()) {
    return;
  }
  int ntraces = traces.size();

  int iter = 0;
  std::vector<Eigen::Vector2d> cur_line_edges;
  std::vector<Eigen::Vector3d> cur_points;
  std::vector<Eigen::Vector3d> cur_colors;

  for (int tid = 0; tid < ntraces; tid++) {
    int nsteps = traces.at(tid).stream_pts_.size();

    for (int i = 0; i < nsteps-1; i++ ) {
      Eigen::Vector3d edge = traces.at(tid).stream_pts_[i].start_pt_ - traces.at(tid).stream_pts_[i + 1].start_pt_;

      if(edge[2] < 0) {
        edge *= -1;
      }

      Eigen::Vector3d rgb_color = SegmentColor(edge, rot_mat);
      cur_points.push_back( traces.at(tid).stream_pts_[i].start_pt_ );
      cur_colors.push_back(rgb_color);
    }
    cur_points.push_back( traces.at(tid).stream_pts_[nsteps - 1].start_pt_ );

    for (int i = 0; i < nsteps - 1; i++ ) {
      cur_line_edges.push_back(Eigen::Vector2d(iter+i, iter+i+1) );
    }
    iter = iter + nsteps;
  }

  polyscope::CurveNetwork* streamlines;

  streamlines = polyscope::registerCurveNetwork(name + "streamline", cur_points, cur_line_edges);
  streamlines->addEdgeColorQuantity("direction colors", cur_colors);
  streamlines->setTransparency(1);
  streamlines->setRadius(0.002);
}

static Eigen::Matrix3d Euler2RotMat(const double roll,
                                    const double pitch,
                                    const double yaw )
{
  Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

  Eigen::Quaterniond q = yawAngle * pitchAngle * rollAngle;
  return q.matrix();
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
int sample_density = 100;
double stream_pt_eps = 1e-2;

int iso_surface_idx = -1;
Eigen::Vector3d rot_axis(0, 0, 1);
double rot_angle = 0;

Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();

ParametrizationType param_type = kSeamless;
StreamLineType streamline_type = kInitPerp;
StreamLineTracingType streamline_tracing_type = kRandomCentroid;

std::string save_folder;

std::vector<CubeCover::Streamline> input_traces, grad_traces, dual_traces, best_match_traces;

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

  if(ImGui::CollapsingHeader("Streamlines tracing", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Combo("Streamline Type", (int*)&streamline_type, "Init\0Gradient\0Init Perp\0Best Match Init\0");
    ImGui::Combo("Tracing Pt Type", (int*)&streamline_tracing_type, "Random Sample Centroid\0Random\0Grid Points\0");
    ImGui::InputInt("Sample Density", &sample_density);
    ImGui::InputDouble("Stream Point Eps", &stream_pt_eps);
    if(ImGui::Button("Trace Stream lines")) {
      std::vector<CubeCover::Streamline> traces;
      int max_iter_per_trace = 700;

      Eigen::MatrixXd frame_vecs;
      std::string frame_name = "input ";
      std::vector<Eigen::MatrixXd> frame_list;

      auto pc_mesh = polyscope::getPointCloud("centroid pc");

      switch (streamline_type) {
      case kInit: {
        frame_vecs = frames;
        frame_list = ExtractFrameVectors(T.rows(), frame_vecs);
        if (streamline_tracing_type == kGridPt) {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces, stream_pt_eps);
        }
        else {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density, streamline_tracing_type == kRandom, stream_pt_eps);
        }
        
        RenderStreamlines(traces, V, T, 1, frame_name, rot_mat);
        input_traces = std::move(traces);
        break;
      }
      case kGradient: {
        frame_vecs = ComputeGradient(V, mesh, values);
        frame_name = "gradient ";
        frame_list = ExtractFrameVectors(T.rows(), frame_vecs);
        if (streamline_tracing_type == kGridPt) {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces, stream_pt_eps);
        }
        else {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density, streamline_tracing_type == kRandom, stream_pt_eps);
        }
        RenderStreamlines(traces, V, T, 1, frame_name, rot_mat);
        grad_traces = std::move(traces);
        break;
      }
      case kInitPerp:{
        frame_vecs = frames_to_trace;
        frame_name = "dual input ";
        frame_list = ExtractFrameVectors(T.rows(), frame_vecs);
        if (streamline_tracing_type == kGridPt) {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces, stream_pt_eps);
        }
        else {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density, streamline_tracing_type == kRandom, stream_pt_eps);
        }
        RenderStreamlines(traces, V, T, 1, frame_name, rot_mat);
        dual_traces = std::move(traces);
        break;
      }
      case kInitBestMatchGrad: {
        Eigen::MatrixXd grad_vec = ComputeGradient(V, mesh, values);
        std::vector<Eigen::MatrixXd> grad_list = ExtractFrameVectors(T.rows(), grad_vec);
        frame_vecs = frames_to_trace;
        frame_list = GetBestMatchFrames(grad_list, frames);
        frame_name = "best match input ";
        if (streamline_tracing_type == kGridPt) {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, values, max_iter_per_trace, traces, stream_pt_eps);
        }
        else {
            CubeCover::TraceStreamlines(V, mesh, frame_vecs, max_iter_per_trace, traces, sample_density, streamline_tracing_type == kRandom, stream_pt_eps);
        }
        RenderStreamlines(traces, V, T, 1, frame_name, rot_mat);
        best_match_traces = std::move(traces);
      }
      }

      for(int i = 0; i < frame_list.size(); i++) {
        pc_mesh->addVectorQuantity(frame_name + std::to_string(i), frame_list[i]);
      }

    }
  }

  if (ImGui::CollapsingHeader("Visualization Options", ImGuiTreeNodeFlags_DefaultOpen)) {
    if(ImGui::Button("Iso Lines extraction")) {
      Eigen::MatrixXd P;
      Eigen::MatrixXi E;

      Eigen::MatrixXd P2;
      Eigen::MatrixXi E2;

      CubeCover::ExtractIsolines(V, mesh, values, P, E, P2, E2);

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

    if(ImGui::Button("Compute Gradient")) {
      Eigen::MatrixXd grad = ComputeGradient(V, mesh, values);
      std::vector<Eigen::MatrixXd> grad_vec = ExtractFrameVectors(T.rows(), grad);

      std::vector<Eigen::VectorXd> errs = GetFrameDifference(frames, grad, grad_vec.size());

      auto tet_mesh = polyscope::getVolumeMesh("tet soup mesh");
      auto pc_mesh = polyscope::getPointCloud("centroid pc");
      for(int i = 0; i < 3; i++) {
        pc_mesh->addVectorQuantity("gradient " + std::to_string(i), grad_vec[i]);
        tet_mesh->addCellScalarQuantity("error " + std::to_string(i), errs[i]);
      }
    }

  }

  if (ImGui::CollapsingHeader("Euler Angles", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::InputDouble("axis 0", &rot_axis[0]);
    ImGui::InputDouble("axis 1", &rot_axis[1]);
    ImGui::InputDouble("axis 2", &rot_axis[2]);
    ImGui::InputDouble("rot angle", (&rot_angle));

    if(ImGui::Button("Compute Rot Mat and update render")) {
      std::cout << rot_angle << std::endl;
      rot_mat = Eigen::AngleAxis<double>(rot_angle / 360.0 * 2 * M_PI, rot_axis.normalized());
      std::cout << rot_mat << std::endl;

      RenderStreamlines(input_traces, V, T, 1, "input ", rot_mat);
      RenderStreamlines(grad_traces, V, T, 1, "gradient ", rot_mat);
      RenderStreamlines(dual_traces, V, T, 1, "dual input ", rot_mat);
      RenderStreamlines(best_match_traces, V, T, 1, "best match input ", rot_mat);

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
    std::cerr << "Usage: isosurfaceviewer (.mesh file) (.bfra file or .fra file) [.hexex file]" << std::endl;
    return -1;
  }

  std::string meshfile = argv[1];
  std::string frafile = argv[2];
  std::string permfile;

  std::string hexexfile = "";
  if (argc == 4) {
      hexexfile = argv[3];
  }
      

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

  if (hexexfile == "") {
      IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling);
  }
  else {
      Eigen::MatrixXd V1;
      Eigen::MatrixXi T1;
      if (!CubeCover::readHexEx(hexexfile, V1, T1, values))
      {
          std::cerr << "error reading the .hexex file" << std::endl;
          return -1;
      }
      else {
          if (V1.rows() != V.rows() || (T - T1).norm() != 0) {
              std::cout << "mesh dismatched!" << std::endl;
              return -1;
          }
      }
  }

  

  std::vector<Eigen::MatrixXd> frame_vec = ExtractFrameVectors(T.rows(), frames);
  stream_pt_eps = 0;

  polyscope::init();

  auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
  psMesh->setTransparency(0.2);
  psMesh->setEnabled(false);

  Eigen::MatrixXd soup_V;
  Eigen::MatrixXi soup_T;

  std::tie(soup_V, soup_T) = GetTetSoup(V, T);

  auto tet_mesh = polyscope::registerTetMesh("tet soup mesh", soup_V, soup_T);

  std::vector<Eigen::Vector3d> centroids;
  for(int i = 0; i < T.rows(); i++) {
    Eigen::Vector3d c;
    c.setZero();
    for(int j = 0; j < 4; j++) {
      c += V.row(T(i, j));
    }
    c /= 4;
    centroids.emplace_back(c);
  }

  auto pc_mesh = polyscope::registerPointCloud("centroid pc", centroids);
  pc_mesh->setEnabled(false);

  std::vector<int> face_ids = {};
  for(int i = 0; i < 3; i++) {
    pc_mesh->addVectorQuantity("frame " + std::to_string(i), frame_vec[i]);
  }

  RenderScalarFields(tet_mesh, values);

  // Add the callback
  polyscope::state::userCallback = callback;

  // visualize!
  polyscope::show();
}
