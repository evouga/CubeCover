#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include "ReadFrameField.h"
#include <iostream>
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include "SingularCurveNetwork.h"
#include "readMeshFixed.h"
#include "polyscope/surface_mesh.h"

int main(int argc, char *argv[])
{
    if (argc != 3 && argc != 4)
    {
        std::cerr << "Usage: singularityviewer (.mesh file) (.fra file) [.perm file]" << std::endl;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    
    std::string meshfile = argv[1];
    std::string frafile = argv[2];
    std::string permfile;
    if (argc == 4)
        permfile = argv[3];

    Eigen::MatrixXi F;
    if (!CubeCover::readMESH(meshfile, V, T, F))
        return -1;

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    if (!CubeCover::readFrameField(frafile, permfile, T, frames, assignments, true))
        return -1;

    CubeCover::TetMeshConnectivity mesh(T);

    CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return -1;

    Eigen::MatrixXd P;
    Eigen::MatrixXi E;
    CubeCover::extractSingularCurveNetwork(V, mesh, *field, P, E);

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

    polyscope::registerCurveNetwork("Singular Curves", P, E);
    auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
    psMesh->setTransparency(0.2);

    // visualize!
    polyscope::show();

    delete field;
}