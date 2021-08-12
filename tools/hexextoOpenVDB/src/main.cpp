#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "VoxelUtils.h"
#include "SceneInfo.h"
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>


/*

This file contains a basic export tool to VDB files 
 - with some amount of scriptable configuration for controlling the final 3d texture. 
 - exports line emmitters and fog.     


*/



int main(int argc, char *argv[])
{
 /*   if (argc > 3 || argc < 2 )
    {
        std::cerr << "Usage: hexextoOpenVDB (.hexex file) [bad_verts path]" << std::endl;
        return -1;
    }
*/
    bool showViz = true;
    std::string badverts;

    if (argc >= 3)
    {
        badverts = argv[2];
        if (badverts != "1")
            showViz = false;
    }


// Config settings.  
    double cells = 2.;
    // double cell_res = 64.;
    double cell_res = 32.;
    double sample_res = cells * cell_res;  // target_cells * res_per_cell 
    double line_w = .05;
    double border_w = .2; // 
    // double cosmic_background = 0.0000; // this is the glow of the parameterisation.   Try setting to like .04
    embedding cur_embed = embedding::PARAM_SPACE;  // PARAM_SPACE



    if( cur_embed == embedding::PARAM_SPACE )
    {
        cell_res = cell_res / 8.;


    }





//    std::string hexexfile = argv[1];
    // std::string hexexfile = "~/Documents/MATLAB/integrable-frames-3d/output_frames_dir/notch5_500.hexex";
     // std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/notch5_500.hexex";
  //  std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/sphere_r0.17.hexex";
    std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/triangular_bipyramid_int.hexex";
    std::string permfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/triangular_bipyramid.perm";

    SceneInfo sc(hexexfile, sample_res);

    // just in case for debugging.
   /* std::ifstream  src(hexexfile, std::ios::binary);
    std::ofstream  dst("./prevToVol.hexex",   std::ios::binary);

    dst << src.rdbuf();
*/




	// draw tet mesh

    double BIG_NUM = 100000000000000000.0;

	int ntets = sc.T.rows();
    int nverts = sc.V.rows();
	Eigen::MatrixXd P_viztets(4 * ntets, 3);
	Eigen::MatrixXi E_viztets(6 * ntets, 2);



 



    
    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid_r = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_g = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_b = openvdb::FloatGrid::create();



    openvdb::FloatGrid::Ptr grid_smoke_r = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_smoke_g = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_smoke_b = openvdb::FloatGrid::create();



    openvdb::FloatGrid::Ptr grid_strength = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_smoke_density = openvdb::FloatGrid::create();
    openvdb::Vec3SGrid::Ptr grid_smoke_color = openvdb::Vec3SGrid::create();
    // openvdb::Vec3SGrid::Ptr grid_smoke_color2 = openvdb::Vec3SGrid::create();    


    grid_r->setName("emission_r");
    grid_g->setName("emission_g");
    grid_b->setName("emission_b");
    grid_smoke_r->setName("smoke_r");
    grid_smoke_g->setName("smoke_g");
    grid_smoke_b->setName("smoke_b");

    grid_strength->setName("emission_strength");
    grid_smoke_density->setName("smoke_density");
    grid_smoke_color->setName("smoke_color");
    // grid_smoke_color2->setName("smoke_color2");

    grid_r->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_g->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_b->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_r->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_g->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_b->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_strength->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_density->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_color->setGridClass(openvdb::GRID_FOG_VOLUME);
    // grid_smoke_color2->setGridClass(openvdb::GRID_FOG_VOLUME);

    openvdb::Coord ijk;

    // Draw voxels
    
    // stampLatticeView(sc,
    //                  acc_r, 
    //                  acc_g, 
    //                  acc_b, 
    //                  acc_smoke_r, 
    //                  acc_smoke_g, 
    //                  acc_smoke_b, 
    //                  acc_strength, // for now blackbody and emission are the same
    //                  acc_smoke_density,
    //                  acc_smoke_color);


    // sc.V_curr = embedding::WORLD_SPACE;
    stampParamView(sc);
    


// Show vertex positions
    // for (int i = 0; i < nverts; i++)
    // {
    //     Eigen::Vector3d round = V_pixel.row(i);
    //     openvdb::Coord ijk;
    //     ijk[0] = int( round(0) );
    //     ijk[1] = int( round(1) );
    //     ijk[2] = int( round(2) );
    //     acc_strength.setValue(ijk, float( 1. )); 
    //     acc_r.setValue(ijk, 1.);
    //     acc_g.setValue(ijk, 1.);
    //     acc_b.setValue(ijk, 1.);

    // }

    grid_r->setName("emission_r");
    grid_g->setName("emission_g");
    grid_b->setName("emission_b");
    grid_smoke_r->setName("smoke_r");
    grid_smoke_g->setName("smoke_g");
    grid_smoke_b->setName("smoke_b");
    grid_strength->setName("emission_strength");
    grid_smoke_density->setName("smoke_density");
    grid_smoke_color->setName("smoke_color");

    openvdb::GridPtrVecPtr grids(new openvdb::GridPtrVec);
    grids->push_back(grid_r);
    grids->push_back(grid_g);
    grids->push_back(grid_b);
    grids->push_back(grid_smoke_r);
    grids->push_back(grid_smoke_g);
    grids->push_back(grid_smoke_b);
    grids->push_back(grid_strength);
    grids->push_back(grid_smoke_density);
    grids->push_back(grid_smoke_color);

    std::ofstream ofile("mygrids.vdb", std::ios_base::binary);
    openvdb::io::Stream(ofile).write(*grids);

        // Print all active ("on") voxels by means of an iterator.
    // for (openvdb::FloatGrid::ValueOnCIter iter = grid_smoke_density->cbeginValueOn(); iter; ++iter) {
    //     std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
    // }

}

/*
    for ( int i = 0; i < sample_rngbnd(0); i++ )
    {
        for ( int j = 0; j < sample_rngbnd(1); j++ )
        {
            for ( int k = 0; k < sample_rngbnd(2); k++ )
            {
                Eigen:Vector3d samp_point(i, j, k);
                samp_point = sample_rngbnd * samp_point + mesh_min;
            }       
        }
    }


	for (int i = 0; i < ntets; i++)
	{
		int idx = 0;
		for (int j = 0; j < 4; j++)
		{
			P_viztets.row(4 * i + j) = param.row(4 * i + j);
			
			for (int k = j + 1; k < 4; k++)
			{
				E_viztets(6 * i + idx, 0) = 4 * i + j;
				E_viztets(6 * i + idx, 1) = 4 * i + k;
				idx++;
			}
		}

	}
*/

/*

	//

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


    extractIsolines(V, mesh, param, P, E, P2, E2, badvert_ids);
        std::cout << "BAD VERTS SIZE " << badvert_ids.size() <<std::endl;

    extractPoints(V, mesh, param, points, colors);


    for (size_t i = 0; i < V.rows(); i++) {
      tet_centers.push_back( glm::vec3{V(i,0), V(i,1), V(i,2)} );
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
*/

/*

        // Following this idea: http://steve.hollasch.net/cgindex/geometry/ptintet.html
        // if this is slow, can move to python and call this vectorized code 
  Eigen::MatrixXd sgn = MatrixXd::Constant(4,4,1.);
        sgn.block<3,1>(0,0) = A;
        sgn.block<3,1>(1,0) = B;
        sgn.block<3,1>(2,0) = C;
        sgn.block<3,1>(3,0) = D;

        // switch sign so it's positive
        if ( sgn.determinant() < 0. )
        {
            Eigen::Vector3d A = V_pixel.row(T[1]);
            Eigen::Vector3d B = V_pixel.row(T[0]);
        }


*/