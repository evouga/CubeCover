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

    CubeCover::CubeCoverOptions opt;
    opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_INTEGERGRID;
    opt.assignmentHandling = (argc == 5 ? CubeCover::CubeCoverOptions::AssignmentHandling::AH_USEPROVIDED : CubeCover::CubeCoverOptions::AssignmentHandling::AH_RECOMPUTE);    
    opt.boundaryConditions = CubeCover::CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER;
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