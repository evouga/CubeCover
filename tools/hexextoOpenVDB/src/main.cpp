#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "ExtractIsolines.h"
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

#include <fstream>

int main(int argc, char *argv[])
{
 /*   if (argc > 3 || argc < 2 )
    {
        std::cerr << "Usage: hexextoOpenVDB (.hexex file) [bad_verts path]" << std::endl;
        return -1;
    }
*/
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

//    std::string hexexfile = argv[1];
    // std::string hexexfile = "~/Documents/MATLAB/integrable-frames-3d/output_frames_dir/notch5_500.hexex";
    std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/notch5_500.hexex";

    // just in case for debugging.
   /* std::ifstream  src(hexexfile, std::ios::binary);
    std::ofstream  dst("./prevToVol.hexex",   std::ios::binary);

    dst << src.rdbuf();
*/
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
    bool showViz = true;
    std::string badverts;

    if (argc >= 3)
    {
        badverts = argv[2];
        if (badverts != "1")
            showViz = false;
    }

    std::cout << badverts << std::endl;



    CubeCover::TetMeshConnectivity mesh(T);
    
    Eigen::MatrixXd P;
    Eigen::MatrixXi E;

    Eigen::MatrixXd P2;
    Eigen::MatrixXi E2;

    std::vector<glm::vec3> points;
    std::vector<std::array<double, 3>> colors;

    std::vector<glm::vec3> tet_centers;
    std::vector<int> badvert_ids;


    extractIsolines(V, mesh, values, P, E, P2, E2, badvert_ids);
        std::cout << "BAD VERTS SIZE " << badvert_ids.size() <<std::endl;

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

    if (showViz)
    {

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
    else
    {
        std::ofstream ofs(badverts);
        if (!ofs)
        {
            return false;
        }
        int nbadverts = badvert_ids.size();
        // ofs << "ids " << V.rows() << " " << nbadverts << std::endl;
        // for (int i = 0; i < nbadverts; i++)
        // {
        //     ofs << badvert_ids[i] << std::endl;
        // }
        std::cout << T.rows() << nbadverts<< std::endl;
        ofs << "ids " << T.rows() << " " << nbadverts << std::endl;
        for (int i = 0; i < nbadverts; i++)
        {
            ofs << badvert_ids[i] << std::endl;
        }
    }

}
