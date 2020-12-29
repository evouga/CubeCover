#include "readMeshFixed.h"
#include <iostream>
#include "igl/writeOBJ.h"
#include "SurfaceExtraction.h"
#include "ReadFrameField.h"
#include "CubeCover.h"
#include "TetMeshConnectivity.h"

int main()
{    
    //std::string example = "three";
    std::string example = "tetrahedron_100";
    
    std::string meshfile = "meshes/" + example + ".mesh";
    std::string framefile = "meshes/" + example + ".fra";
    std::string permfile = "meshes/" + example + ".perm";

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXi F;
    if (!CubeCover::readMESH("../" + meshfile, V, T, F))
    {
        if (!CubeCover::readMESH("../../" + meshfile, V, T, F))
        {
            std::cerr << "Cannot open mesh file" << std::endl;
            return -1;
        }
    }
    CubeCover::TetMeshConnectivity mesh(T);

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    
    //if (!CubeCover::readFrameField("../" + framefile, "../" + permfile, T, frames, assignments, true))
    if (!CubeCover::readFrameField("../" + framefile, "", T, frames, assignments, true))
    {        
        //if (!CubeCover::readFrameField("../../" + framefile, "../../" + permfile, T, frames, assignments, true))
        if (!CubeCover::readFrameField("../../" + framefile, "", T, frames, assignments, true))
        {
            std::cerr << "Cannot open frame field file" << std::endl;
            return -1;
        }
    }

    Eigen::MatrixXd soupValues;
    CubeCover::CubeCoverOptions opt;
    opt.parameterizationType = CubeCover::CubeCoverOptions::ParameterizationType::PT_SEAMLESS;
    opt.verbose = true;
    CubeCover::cubeCover(V, T, frames, assignments, soupValues, opt);
    
    Eigen::MatrixXd isoV;
    Eigen::MatrixXi isoF;
    CubeCover::isosurfaceSoup(V, mesh, soupValues, isoV, isoF);

    igl::writeOBJ("iso.obj", isoV, isoF);
}