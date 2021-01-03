#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "ExtractIsolines.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: isolineviewer (.hexex file)" << std::endl;
        return -1;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    std::string hexexfile = argv[1];
    
    Eigen::MatrixXd values;
    if (!CubeCover::readHexEx(hexexfile, V, T, values))
    {
        std::cerr << "error reading the .hexex file" << std::endl;
        return -1;
    }

    CubeCover::TetMeshConnectivity mesh(T);

    
    Eigen::MatrixXd P;
    Eigen::MatrixXi E;

    extractIsolines(V, mesh, values, P, E);
    
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

    auto *psCurves = polyscope::registerCurveNetwork("Isolines", P, E);
    psCurves->setRadius(0.01);
    auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
    psMesh->setTransparency(0.2);
    

    // visualize!
    polyscope::show();
}
