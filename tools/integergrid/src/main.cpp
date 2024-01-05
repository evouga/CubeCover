#include <iostream>
#include <Eigen/Core>
#include "ReadFrameField.h"
#include "readMeshFixed.h"
#include "CubeCover.h"
#include "WriteHexEx.h"

int main(int argc, char* argv[])
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: integergrid (.mesh file) (.fra file) [.perm file] (output file)" << std::endl;
        return -1;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    std::string meshfile = argv[1];
    std::string frafile = argv[2];
    std::string permfile;
    if (argc == 5)
        permfile = argv[3];

    std::string outfile;
    if (argc == 5)
        outfile = argv[4];
    else
        outfile = argv[3];

    Eigen::MatrixXi F;
    if (!CubeCover::readMESH(meshfile, V, T, F))
    {
        std::cerr << "could not read .mesh file " << meshfile << std::endl;
        return -1;
    }

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
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

    CubeCover::CubeCoverOptions opt;
    opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_INTEGERGRID;
//    opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_SEAMLESS;
    opt.assignmentHandling = (argc == 5 ? CubeCover::CubeCoverOptions::AssignmentHandling::AH_USEPROVIDED : CubeCover::CubeCoverOptions::AssignmentHandling::AH_RECOMPUTE);
    opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER;

//    opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FREE;
    opt.solver = CubeCover::CubeCoverOptions::MIPSolver::MS_COMISO;

    // set to something non-zero if you want curl-correction. 1.0 == 100% change in the input frames allowed.
    opt.curlCorrection = 0.0;

    opt.verbose = true;
    Eigen::MatrixXd values;
    if (!CubeCover::cubeCover(V, T, frames, assignments, values, opt))
    {        
        return -1;
    }

    if (!CubeCover::writeHexEx(outfile, V, T, values))
    {
        std::cerr << "could not write output file " << outfile << std::endl;
        return -1;
    }

    return 0;
}