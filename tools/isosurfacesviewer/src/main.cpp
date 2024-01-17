#include <Eigen/Core>
#include <Eigen/Dense>
#include <filesystem>
#include <iostream>

#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/writeMESH.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"

#include "TetMeshConnectivity.h"
#include "ReadFrameField.h"
#include "readMeshFixed.h"
#include "ReadHexEx.h"
#include "WriteHexEx.h"
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

static Eigen::MatrixXd NormalizePtsGivenCenterAndScalingRatio(const Eigen::MatrixXd& pos, Eigen::RowVector3d& center, double& scaling_ratio) {
  int npts = pos.rows();
  Eigen::MatrixXd normalize_pts = pos;
  for(int i = 0; i < npts; i++) {
    normalize_pts.row(i) = (normalize_pts.row(i) - center) / scaling_ratio;
  }
  return normalize_pts;
}

// normalize point clouds to [-1, 1] x [-1, 1] x [-1, 1]
static Eigen::MatrixXd NormalizePts(const Eigen::MatrixXd& pos, Eigen::RowVector3d& center, double& scaling_ratio) {
  Eigen::RowVector3d min_corner, max_corner;
  min_corner = pos.colwise().minCoeff().transpose();
  max_corner = pos.colwise().maxCoeff().transpose();

  center = (min_corner + max_corner) / 2;
  scaling_ratio = (max_corner - min_corner).maxCoeff();

  return NormalizePtsGivenCenterAndScalingRatio(pos, center, scaling_ratio);

}

static std::pair<Eigen::MatrixXd, Eigen::MatrixXi> ExportBoundaryMesh(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh) {
	// make a mesh out of all of the boundary faces
	int nbdry = 0;
	int nfaces = mesh.nFaces();
	for (int i = 0; i < nfaces; i++) {
		if (mesh.isBoundaryFace(i))
			nbdry++;
	}
	Eigen::MatrixXd bdry_V(V.rows(), 3);
	Eigen::MatrixXi bdry_F(nbdry, 3);

	std::unordered_map<int, int> vid_to_bryid;

	int curidx = 0;
	int curvidx = 0;
	for (int i = 0; i < nfaces; i++) {
		if (mesh.isBoundaryFace(i)) {
			for (int j = 0; j < 3; j++) {
				int vid = mesh.faceVertex(i, j);
				if (vid_to_bryid.count(vid)) {
					bdry_F(curidx, j) = vid_to_bryid[vid];
				}
				else {
					vid_to_bryid[vid] = curvidx;
					bdry_F(curidx, j) = curvidx;
					bdry_V.row(curvidx++) = V.row(vid);
				}
			}
			// fix triangle orientations
			int tet = mesh.faceTet(i, 1);
			if (tet == -1) {
				std::swap(bdry_F(curidx, 0), bdry_F(curidx, 1));
			}
			curidx++;
		}
	}
	bdry_V.conservativeResize(curvidx, Eigen::NoChange);
	return { bdry_V, bdry_F };
}

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
  bool is_success = CubeCover::cubeCover(V, T, to_round, assignments, values, opt);

//  int ntets = T.rows();
//  CubeCover::TetMeshConnectivity mesh(T);
//  for(int tid = 0; tid < ntets; tid++) {
//        std::array<int, 3> bnd_vert_id_in_tet = {-1, -1, -1};
//        for(int i = 0; i < 4; i++) {
//          int fid = mesh.tetFace(tid, i);
//          if(mesh.isBoundaryFace(fid)) {
//                for(int k = 0; k < 3; k++) {
//                  int vid = mesh.faceVertex(fid, k);
//                  for(int vi = 0; vi < 4; vi++) {
//                    if(mesh.tetVertex(tid, vi) == vid) {
//                      for(int c = 0; c < 3; c++) {
//                        int val = std::round(values(4 * tid + vi, c));
//                        if(std::abs(values(4 * tid + vi, c) - val) < 1e-3) {
//                          values(4 * tid + vi, c) = val;
//                        }
//                      }
//                    }
//                  }
//                }
//                break;
//          }
//        }
//  }

  return is_success;
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

// Helper functions
double cc(double na, double nd) {
  return na * std::cos(nd * M_PI / 180.0);
}

double ss(double na, double nd) {
  return na * std::sin(nd * M_PI / 180.0);
}

// reference: https://github.com/fury-gl/fury/blob/2a35d023138b99cb336247d3e893780f9a55e6ac/fury/colormap.py#L62
Eigen::Vector3d Boys2RBG(const Eigen::Vector3d& v) {
  Eigen::Vector3d color;

  Eigen::Vector3d norm = v.normalized();

  double x = norm[0], y = norm[1], z = norm[2];

  // Calculations similar to Python code
  double x2 = x * x, y2 = y * y, z2 = z * z;
  double x3 = x * x2, y3 = y * y2, z3 = z * z2;
  double z4 = z2 * z2;
  double xy = x * y, xz = x * z, yz = y * z;

  // Constants
  double hh1 = 0.5 * (3 * z2 - 1) / 1.58;
  double hh2 = 3 * xz / 2.745;
  double hh3 = 3 * yz / 2.745;
  double hh4 = 1.5 * (x2 - y2) / 2.745;
  double hh5 = 6 * xy / 5.5;
  double hh6 = (1 / 1.176) * 0.125 * (35 * z4 - 30 * z2 + 3);
  double hh7 = 2.5 * x * (7 * z3 - 3 * z) / 3.737;
  double hh8 = 2.5 * y * (7 * z3 - 3 * z) / 3.737;
  double hh9 = ((x2 - y2) * 7.5 * (7 * z2 - 1)) / 15.85;
  double hh10 = ((2 * xy) * (7.5 * (7 * z2 - 1))) / 15.85;
  double hh11 = 105 * (4 * x3 * z - 3 * xz * (1 - z2)) / 59.32;
  double hh12 = 105 * (-4 * y3 * z + 3 * yz * (1 - z2)) / 59.32;

  double s0 = -23.0;
  double s1 = 227.9;
  double s2 = 251.0;
  double s3 = 125.0;

  double ss23 = ss(2.71, s0);
  double cc23 = cc(2.71, s0);
  double ss45 = ss(2.12, s1);
  double cc45 = cc(2.12, s1);
  double ss67 = ss(0.972, s2);
  double cc67 = cc(0.972, s2);
  double ss89 = ss(0.868, s3);
  double cc89 = cc(0.868, s3);

  double X = 0.0;

  X = X + hh2 * cc23;
  X = X + hh3 * ss23;

  X = X + hh5 * cc45;
  X = X + hh4 * ss45;

  X = X + hh7 * cc67;
  X = X + hh8 * ss67;

  X = X + hh10 * cc89;
  X = X + hh9 * ss89;

  double Y = 0.0;

  Y = Y + hh2 * -ss23;
  Y = Y + hh3 * cc23;

  Y = Y + hh5 * -ss45;
  Y = Y + hh4 * cc45;

  Y = Y + hh7 * -ss67;
  Y = Y + hh8 * cc67;

  Y = Y + hh10 * -ss89;
  Y = Y + hh9 * cc89;

  double Z = 0.0;

  Z = Z + hh1 * -2.8;
  Z = Z + hh6 * -0.5;
  Z = Z + hh11 * 0.3;
  Z = Z + hh12 * -2.5;

  double w_x = 4.1925, trl_x = -2.0425;
  double w_y = 4.0217, trl_y = -1.8541;
  double w_z = 4.0694, trl_z = -2.1899;

  color[0] = 0.9 * std::abs(((X - trl_x) / w_x)) + 0.05;
  color[1] = 0.9 * std::abs(((Y - trl_y) / w_y)) + 0.05;
  color[2] = 0.9 * std::abs(((Z - trl_z) / w_z)) + 0.05;

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

  return Boys2RBG(dir_hat);
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

static void GetStreamlines(const std::vector<CubeCover::Streamline>& traces, Eigen::MatrixXd& p_start, Eigen::MatrixXd& p_end, Eigen::MatrixXd& colors, Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity()) {
	if (traces.empty()) {
		return;
	}
	int ntraces = traces.size();

	p_start.resize(10000, 3);
	p_end.resize(10000, 3);
	colors.resize(10000, 4);

	int id = 0;

	for (int tid = 0; tid < ntraces; tid++) {
		int nsteps = traces.at(tid).stream_pts_.size();

		for (int i = 0; i < nsteps - 1; i++) {
			Eigen::Vector3d edge = traces.at(tid).stream_pts_[i].start_pt_ - traces.at(tid).stream_pts_[i + 1].start_pt_;
			Eigen::Vector3d rgb_color = SegmentColor(edge, rot_mat);

			if (id >= p_start.rows()) {
				p_start.conservativeResize(p_start.rows() + 10000, Eigen::NoChange);
				p_end.conservativeResize(p_end.rows() + 10000, Eigen::NoChange);
				colors.conservativeResize(colors.rows() + 10000, Eigen::NoChange);
			}

			p_start.row(id) = traces[tid].stream_pts_[i].start_pt_.transpose();
			p_end.row(id) = traces[tid].stream_pts_[i + 1].start_pt_.transpose();
			colors.row(id) << rgb_color[0], rgb_color[1], rgb_color[2], 1;
			id++;
		}
	}
	p_start.conservativeResize(id, Eigen::NoChange);
	p_end.conservativeResize(id, Eigen::NoChange);
	colors.conservativeResize(id, Eigen::NoChange);
}

static void SaveStreamlines(const std::string file, const Eigen::MatrixXd& p_start, const Eigen::MatrixXd& p_end, const Eigen::MatrixXd& colors) {
	std::ofstream ofs(file);
	for (int i = 0; i < p_start.rows(); i++) {
		ofs << p_start.row(i) << " " << p_end.row(i) << " " << colors.row(i) << std::endl;
	}
}



static Eigen::Matrix3d Euler2RotMat(const double roll, const double pitch, const double yaw )
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

double global_rescaling = 10;
int sample_density = 100;
double stream_pt_eps = 1e-2;

int iso_surface_idx = -1;
Eigen::Vector3d rot_axis(0, 0, 1);
double rot_angle = 0;

Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();

ParametrizationType param_type = kIntegerGrid;
StreamLineType streamline_type = kInitPerp;
StreamLineTracingType streamline_tracing_type = kRandomCentroid;

std::string save_folder;

std::vector<CubeCover::Streamline> input_traces, grad_traces, dual_traces, best_match_traces;


void SaveforRender(std::string folder) {
	if (!std::filesystem::exists(folder)) {
		std::filesystem::create_directory(folder);
	}
	Eigen::MatrixXd bnd_V;
	Eigen::MatrixXi bnd_F;
	std::tie(bnd_V, bnd_F) = ExportBoundaryMesh(V, mesh);

        double scaling_ratio;
        Eigen::RowVector3d center;
        bnd_V = NormalizePts(bnd_V, center, scaling_ratio);

	igl::writeOBJ(folder + "/boundary_mesh.obj", bnd_V, bnd_F);

	Eigen::MatrixXd p_start, p_end, colors;

	// stream lines
	CubeCover::TraceStreamlines(V, mesh, frames_to_trace, values, 700, dual_traces, stream_pt_eps);
	GetStreamlines(dual_traces, p_start, p_end, colors);

        p_start = NormalizePtsGivenCenterAndScalingRatio(p_start, center, scaling_ratio);
        p_end = NormalizePtsGivenCenterAndScalingRatio(p_end, center, scaling_ratio);
	SaveStreamlines(folder + "/stream_lines.txt", p_start, p_end, colors);

	// isolines
	Eigen::MatrixXd P;
	Eigen::MatrixXi E;

	Eigen::MatrixXd P2;
	Eigen::MatrixXi E2;
	CubeCover::ExtractIsolines(V, mesh, values, P, E, P2, E2);

	Eigen::MatrixXd iso_start_pt((E.rows() + E2.rows()), 3), iso_end_pt((E.rows() + E2.rows()), 3), iso_colors((E.rows() + E2.rows()), 4);

	for (int i = 0; i < E.rows(); i++) {
		int vid0 = E(i, 0);
		iso_start_pt.row(i) = P.row(vid0);
		int vid1 = E(i, 1);
		iso_end_pt.row(i) = P.row(vid1);
		iso_colors.row(i) << 0, 0, 0, 1;
	}

	for (int i = 0; i < E2.rows(); i++) {
		int vid0 = E2(i, 0);
		iso_start_pt.row(E.rows() + i) = P2.row(vid0);
		int vid1 = E2(i, 1);
		iso_end_pt.row(E.rows() + i) = P2.row(vid1);
		iso_colors.row(E.rows() + i) << 0, 0, 0, 1;
	}

        iso_start_pt = NormalizePtsGivenCenterAndScalingRatio(iso_start_pt, center, scaling_ratio);
        iso_end_pt = NormalizePtsGivenCenterAndScalingRatio(iso_end_pt, center, scaling_ratio);

	SaveStreamlines(folder + "/isolines.txt", iso_start_pt, iso_end_pt, iso_colors);
}

void callback() {
  ImGui::PushItemWidth(100);
  if (ImGui::CollapsingHeader("Save Options")) {
	  if (ImGui::Button("Save Hexex")) {
		  std::string file = igl::file_dialog_save();
		  CubeCover::writeHexEx(file, V, T, values);
	  }
	  if (ImGui::Button("Save Boundary Mesh")) {
		  std::string file = igl::file_dialog_save();
		  Eigen::MatrixXd bdryV;
		  Eigen::MatrixXi bdryF;
		  std::tie(bdryV, bdryF) = ExportBoundaryMesh(V, mesh);
		  igl::writeOBJ(file, bdryV, bdryF);
	  }
	  if (ImGui::Button("Save For Render")) {
		  std::string folder = igl::file_dialog_save();
		  SaveforRender(folder);
	  }
  }

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

	if(ImGui::Button("Grid Points")) {
	  std::vector<Eigen::Vector3d> pts;
	  for(int i = 0; i < T.rows(); i++) {
		std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> tmp_res = CubeCover::ComputeGridPts(V, mesh, values, i);
		for(auto& pt : tmp_res) {
		  pts.push_back(pt.second);
		}
	  }
	  polyscope::registerPointCloud("grid points", pts);
	}

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
		frame_vecs /= global_rescaling;
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
		grad_vec /= global_rescaling;
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
	  std::vector<Eigen::Vector3d> isoline_colors;
	  for (int i = 0; i < E.rows(); i++) {
		  Eigen::Vector3d dir = P.row(E(i, 0)) - P.row(E(i, 1));
		  isoline_colors.push_back(SegmentColor(dir));
	  }
	  psCurves->addEdgeColorQuantity("direction colors", isoline_colors);
	  psCurves->setRadius(0.002);


	  isoline_colors.clear();
	  auto *psCurves2 = polyscope::registerCurveNetwork("Bad Isolines", P2, E2);
	  for (int i = 0; i < E2.rows(); i++) {
		  Eigen::Vector3d dir = P2.row(E2(i, 0)) - P2.row(E2(i, 1));
		  isoline_colors.push_back(SegmentColor(dir));
	  }
	  psCurves2->addEdgeColorQuantity("direction colors", isoline_colors);
	  psCurves2->setRadius(0.002);
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
	  grad /= global_rescaling;
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

  if (ImGui::CollapsingHeader("Rotation Matrix", ImGuiTreeNodeFlags_DefaultOpen)) {
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


int main(int argc, char *argv[]) {
  if (argc != 3 && argc != 4) {
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
  if (!CubeCover::readMESH(meshfile, V, T, F)) {
	std::cerr << "could not read .mesh file " << meshfile << std::endl;
	return -1;
  }

  if (!CubeCover::readFrameField(frafile, permfile, T, frames, assignments, true)) {
	std::cerr << "could not read frames/permutations" << std::endl;
	return -1;
  }

  std::cout << "tet mesh: " << T.rows() << ", frame size: " << frames.rows() << std::endl;

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

  Eigen::MatrixXd bdryV;
  Eigen::MatrixXi bdryF;
  std::tie(bdryV, bdryF) = ExportBoundaryMesh(V, mesh);

  if (hexexfile == "") {
	  IntegrateFrames(frames, V, T, param_type, assignments, values, global_rescaling);
  }
  else {
	  Eigen::MatrixXd V1;
	  Eigen::MatrixXi T1;
	  if (!CubeCover::readHexEx(hexexfile, V1, T1, values)) {
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
  stream_pt_eps = 0.2;

  polyscope::init();

  auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", bdryV, bdryF);
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
  polyscope::options::autocenterStructures = false;

  // visualize!
  polyscope::show();
}
