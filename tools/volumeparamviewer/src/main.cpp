#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "ExtractIsolines.h"
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

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


	// draw tet mesh

	int ntets = T.rows();
	Eigen::MatrixXd P_viztets(4 * ntets, 3);
	Eigen::MatrixXi E_viztets(6 * ntets, 2);

	for (int i = 0; i < ntets; i++)
	{
		int idx = 0;
		for (int j = 0; j < 4; j++)
		{
			P_viztets.row(4 * i + j) = values.row(4 * i + j);
			
			for (int k = j + 1; k < 4; k++)
			{
				E_viztets(6 * i + idx, 0) = 4 * i + j;
				E_viztets(6 * i + idx, 1) = 4 * i + k;
				idx++;
			}
		}

	}

	//




    CubeCover::TetMeshConnectivity mesh(T);
    
    Eigen::MatrixXd P;
    Eigen::MatrixXi E;

    Eigen::MatrixXd P2;
    Eigen::MatrixXi E2;

    std::vector<glm::vec3> points;
    std::vector<std::array<double, 3>> colors;

    std::vector<glm::vec3> tet_centers;


    extractIsolines(V, mesh, values, P, E, P2, E2);
    extractPoints(V, mesh, values, points, colors);


    for (size_t i = 0; i < V.rows(); i++) {
      tet_centers.push_back( glm::vec3{V(i,0), V(i,1), V(i,2)} );
    }
    
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


    // auto *psCurves = polyscope::registerCurveNetwork("Isolines", P, E);
    // psCurves->setRadius(0.003);
	auto *psVizTets = polyscope::registerCurveNetwork("viztets", P_viztets, E_viztets);
	psVizTets->setRadius(0.003);
	psVizTets->setEnabled(false);


    auto *psCurves2 = polyscope::registerCurveNetwork("Bad Isolines", P2, E2);
    psCurves2->setRadius(0.003);

    auto *psPoints = polyscope::registerPointCloud("Mesh Verts", points);
    psPoints->setPointRadius(0.003);
    psPoints->addColorQuantity("pos_color", colors)->setEnabled(true);
    psPoints->setEnabled(false);

    auto *psPoints2 = polyscope::registerPointCloud("Mesh Interior", tet_centers);
    psPoints2->setPointRadius(0.003);
    psPoints2->addColorQuantity("pos_color", colors)->setEnabled(true);


    auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
    psMesh->setTransparency(0.2);
    

    // visualize!
    polyscope::show();
}
