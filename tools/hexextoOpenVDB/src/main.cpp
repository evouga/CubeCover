#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "VoxelUtils.h"
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>

#include <fstream>

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
    double cells = 4;
    double cell_res = 64.;
    double sample_res = cells * cell_res;  // target_cells * res_per_cell 
    double line_w = .05;
    double border_w = .2; // 
    double cosmic_background = 0.0000; // this is the glow of the parameterisation.   Try setting to like .04





//    std::string hexexfile = argv[1];
    // std::string hexexfile = "~/Documents/MATLAB/integrable-frames-3d/output_frames_dir/notch5_500.hexex";
     // std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/notch5_500.hexex";
  //  std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/sphere_r0.17.hexex";
    std::string hexexfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/triangular_bipyramid_int.hexex";
    std::string permfile = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/triangular_bipyramid.perm";

    sceneInfo sc(hexexfile, sample_res);

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
    openvdb::FloatGrid::Accessor acc_r = grid_r->getAccessor();
    openvdb::FloatGrid::Accessor acc_g = grid_g->getAccessor();
    openvdb::FloatGrid::Accessor acc_b = grid_b->getAccessor();


    openvdb::FloatGrid::Ptr grid_smoke_r = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_smoke_g = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_smoke_b = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor acc_smoke_r = grid_smoke_r->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_g = grid_smoke_g->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_b = grid_smoke_b->getAccessor();


    openvdb::FloatGrid::Ptr grid_strength = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr grid_smoke_density = openvdb::FloatGrid::create();
    openvdb::Vec3SGrid::Ptr grid_smoke_color = openvdb::Vec3SGrid::create();
    openvdb::Vec3SGrid::Ptr grid_smoke_color2 = openvdb::Vec3SGrid::create();    
    openvdb::FloatGrid::Accessor acc_strength = grid_strength->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_density = grid_smoke_density->getAccessor();
    openvdb::Vec3SGrid::Accessor acc_smoke_color = grid_smoke_color->getAccessor();
    openvdb::Vec3SGrid::Accessor acc_smoke_color2 = grid_smoke_color2->getAccessor();

    grid_r->setName("emission_r");
    grid_g->setName("emission_g");
    grid_b->setName("emission_b");
    grid_smoke_r->setName("smoke_r");
    grid_smoke_g->setName("smoke_g");
    grid_smoke_b->setName("smoke_b");

    grid_strength->setName("emission_strength");
    grid_smoke_density->setName("smoke_density");
    grid_smoke_color->setName("smoke_color");
    grid_smoke_color2->setName("smoke_color2");

    grid_r->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_g->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_b->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_r->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_g->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_b->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_strength->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_density->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_color->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid_smoke_color2->setGridClass(openvdb::GRID_FOG_VOLUME);

    openvdb::Coord ijk;

    // Draw voxels
    for (int t = 0; t < ntets; t++)
    {
        // compute tet bbox 
        Eigen::Vector3d curt_min(BIG_NUM,BIG_NUM,BIG_NUM);
        Eigen::Vector3d curt_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);
        for (int v_idx = 0; v_idx < 4; v_idx++)
        {
            int rowId = sc.T(t, v_idx);
            Eigen::Vector3d cur_point = sc.V_pixel.row( rowId );
            for ( int k = 0; k < 3; k++ )
            {
                if( cur_point(k) < curt_min(k) )
                    curt_min(k) = cur_point(k);
                if( cur_point(k) > curt_max(k) )
                    curt_max(k) = cur_point(k);
            }
        }


// return vertex positions in pixel space.  
        Eigen::Vector3d A = sc.V_pixel.row( sc.T(t, 0) );
        Eigen::Vector3d B = sc.V_pixel.row( sc.T(t, 1) );
        Eigen::Vector3d C = sc.V_pixel.row( sc.T(t, 2) );
        Eigen::Vector3d D = sc.V_pixel.row( sc.T(t, 3) );




        // https://stackoverflow.com/questions/25179693/how-to-check-whether-the-point-is-in-the-tetrahedron-or-not
      
    //    std::cout << curt_min.transpose() << std::endl;
    //    std::cout << curt_max.transpose() << std::endl<< std::endl;
        // iterate over bbox
        openvdb::Coord ijk;
        int &i = ijk[0], &j = ijk[1], &k = ijk[2];
        for (i = (int) curt_min[0]; i < (int) curt_max[0]; ++i) {
            for (j = (int) curt_min[1]; j < (int) curt_max[1]; ++j) {
                for (k = (int) curt_min[2]; k < (int) curt_max[2]; ++k) {
                    Eigen::Vector3d p;
                    p << float(i), float(j), float(k);
                    
                    // Check if a voxel is in a tet ( they in same coordinate system now. )
                    Eigen::VectorXd param_val_tet;
                    bool pIsIn = pointInsideT(A, B, C, D, p, param_val_tet);



            //        std::cout << "return " << param_val_tet.transpose() << std::endl;
             //       std::cout << param_val_tet.transpose() << std::endl << std::endl;
                    if ( pIsIn )
                    {
                        // load tet texture param 
                                        //        acc_smoke_density.setValue(ijk, float( param_val_tet.head(3).sum() )); 
                //        acc_strength.setValue(ijk, float( .5 )); 

                        // param vals 
                        // double s = 3.;
                        // Eigen::Vector3d da = param.row(4 * i + 1)/s - param.row(4 * i)/s;
                        // Eigen::Vector3d db = param.row(4 * i + 2)/s - param.row(4 * i)/s;
                        // Eigen::Vector3d dc = param.row(4 * i + 3)/s - param.row(4 * i)/s;
                        // Eigen::Vector3d orig = param.row(4 * i)/s;
                        


                        // Ambient grid 
                        // Eigen::Vector3d da = B - A;
                        // Eigen::Vector3d db = C - A;
                        // Eigen::Vector3d dc = D - A;
                        // Eigen::Vector3d orig = A;

               
                        Eigen::Vector3d da = sc.param.row(4 * t + 1) - sc.param.row(4 * t);
                        Eigen::Vector3d db = sc.param.row(4 * t + 2) - sc.param.row(4 * t);
                        Eigen::Vector3d dc = sc.param.row(4 * t + 3) - sc.param.row(4 * t);
                        Eigen::Vector3d orig = sc.param.row(4 * t);

                        double s0 = param_val_tet(0);
                        double s1 = param_val_tet(1);
                        double s2 = param_val_tet(2);

                        Eigen::Vector3d interp_param_to_pixel = s0 * da + s1 * db + s2 * dc + orig;
                        // interp_param_to_pixel *= 1. / cells;


                        // This plots the param param per tet in the underlying mesh.  

                        acc_smoke_color.setValue(ijk, openvdb::Vec3s(float( param_val_tet(0) ), 
                                                                     float( param_val_tet(1) ), 
                                                                     float( param_val_tet(2) )) );

                   
                      //   std::cout << "interp_param_to_pixel" << interp_param_to_pixel.transpose() << std::endl<< std::endl;

                        // Find texture coordinate of pixel rounded to a unit cube.  
                        double unit_u = interp_param_to_pixel(0) - floor(interp_param_to_pixel(0));
                        double unit_v = interp_param_to_pixel(1) - floor(interp_param_to_pixel(1));
                        double unit_w = interp_param_to_pixel(2) - floor(interp_param_to_pixel(2));

                       //  std::cout << "round " << interp_param_to_pixel.transpose() << std::endl<< std::endl;


                        // acc_strength.setValue(ijk, float( .04 ));
                        // acc_strength.setValue(ijk, float( cosmic_background ));
                        // acc_r.setValue(ijk, unit_u );
                        // acc_g.setValue(ijk, unit_v );
                        // acc_b.setValue(ijk, unit_w );

                        // This plots the param param per point in parameter space 
                        acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( unit_v ) + float( unit_w ), 
                                                                      float( unit_w ) + float( unit_u ), 
                                                                      float( unit_v ) + float( unit_u )) );

                        // acc_smoke_r.setValue(ijk, unit_v + unit_w );
                        // acc_smoke_g.setValue(ijk, unit_u + unit_w  );
                        // acc_smoke_b.setValue(ijk, unit_u + unit_v  );
                        acc_smoke_r.setValue(ijk, 0.5 );
                        acc_smoke_g.setValue(ijk, 0.5  );
                        acc_smoke_b.setValue(ijk, 0.5  );
                        acc_smoke_density.setValue(ijk, 0. );

                        

                        // Draw the glowing grid.

                        double u_dist = std::min(unit_u, fabs(1. - unit_u) );
                        double v_dist = std::min(unit_v, fabs(1. - unit_v) );
                        double w_dist = std::min(unit_w, fabs(1. - unit_w) );


                        bool border_cell = false; 
                        border_cell = (( unit_u < border_w || unit_u > 1. - border_w )  &&     
                                       ( unit_v < border_w || unit_v > 1. - border_w )) || border_cell;
                        border_cell = (( unit_u < border_w || unit_u > 1. - border_w )  &&     
                                       ( unit_w < border_w || unit_w > 1. - border_w )) || border_cell;
                        border_cell = (( unit_v < border_w || unit_v > 1. - border_w )  &&     
                                       ( unit_w < border_w || unit_w > 1. - border_w )) || border_cell;
                        if ( border_cell )
                        {
                            acc_strength.setValue(ijk, 0.);
                            acc_r.setValue(ijk, 0. );
                            acc_g.setValue(ijk, 0. );
                            acc_b.setValue(ijk, 0. );

                            // acc_smoke_r.setValue(ijk, unit_v + unit_w );
                            // acc_smoke_g.setValue(ijk, unit_u + unit_w  );
                            // acc_smoke_b.setValue(ijk, unit_u + unit_v  );

                            // acc_smoke_r.setValue(ijk, 1. );
                            // acc_smoke_g.setValue(ijk, 1.  );
                            // acc_smoke_b.setValue(ijk, 1.  );

                            acc_smoke_density.setValue(ijk, .0);
                        }
                        else{
                            acc_smoke_density.setValue(ijk, .1);

                            acc_smoke_r.setValue(ijk, .9 );
                            acc_smoke_g.setValue(ijk, .9  );
                            acc_smoke_b.setValue(ijk, .9  );

                            acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( 0.f ), 
                                                                     float( 0.f ), 
                                                                     float( 0.f )) );
                        }


                        if ( v_dist < line_w && w_dist < line_w )
                        {
                            acc_r.setValue(ijk, .6);
                            acc_strength.setValue(ijk, unit_u);
     
                            acc_smoke_density.setValue(ijk, 1.);
                            acc_smoke_r.setValue(ijk, 1. );
                            // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                        }
                        else if ( u_dist < line_w && w_dist < line_w )
                        {
                            acc_g.setValue(ijk, .6);

                            acc_strength.setValue(ijk, unit_v);

                            acc_smoke_density.setValue(ijk, 1.);
                            acc_smoke_g.setValue(ijk, 1.  );
                            // acc_strength.setValue(ijk, float( interp_param_to_pixel(1)) );

                        }
                        else if ( u_dist < line_w && v_dist < line_w )
                        {
                            acc_b.setValue(ijk, .6);

                            acc_strength.setValue(ijk, unit_w);

                            acc_smoke_density.setValue(ijk, 1.);
                            acc_smoke_b.setValue(ijk, 1.  );
                            // acc_strength.setValue(ijk, float( interp_param_to_pixel(2)) );
                        }



                        // if ( ( unit_v < line_w || unit_v > 1. - line_w ) &&     
                        //      ( unit_w < line_w || unit_w > 1. - line_w ) )
                        // {
                        //     acc_r.setValue(ijk, .6);
                        //     acc_strength.setValue(ijk, unit_u);
                        //     acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( unit_u ), 
                        //                                              float( 0.f ), 
                        //                                              float( 0.f )) );
                        //     acc_smoke_density.setValue(ijk, .3);
                        //     acc_smoke_r.setValue(ijk, 1. );
                        //     // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                        // }
                        // else if ( ( unit_u < line_w || unit_u > 1. - line_w ) &&     
                        //         (   unit_w < line_w || unit_w > 1. - line_w ) )
                        // {
                        //     acc_g.setValue(ijk, .6);

                        //     acc_strength.setValue(ijk, unit_v);
                        //     acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( 0.f ), 
                        //                                              float( unit_v ), 
                        //                                              float( 0.f )) );
                        //     acc_smoke_density.setValue(ijk, .3);
                        //     acc_smoke_g.setValue(ijk, 1.  );
                        //     // acc_strength.setValue(ijk, float( interp_param_to_pixel(1)) );

                        // }
                        // else if ( ( unit_u < line_w || unit_u > 1. - line_w ) &&     
                        //         (   unit_v < line_w || unit_v > 1. - line_w ) )
                        // {
                        //     acc_b.setValue(ijk, .6);

                        //     acc_strength.setValue(ijk, unit_w);
                        //     acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( 0.f ), 
                        //                                              float( 0.f ), 
                        //                                              float( unit_w )) );
                        //     acc_smoke_density.setValue(ijk, .3);
                        //     acc_smoke_b.setValue(ijk, 1.  );
                        //     // acc_strength.setValue(ijk, float( interp_param_to_pixel(2)) );
                        // }




                        
/*
                        if ( ( unit_u < line_w || fabs(1. - unit_u) < line_w ) && 
                             ( unit_v < line_w || fabs(1. - unit_v) < line_w ) &&
                             ( unit_w < line_w || fabs(1. - unit_w) < line_w ) )
                        {
                            acc_strength.setValue(ijk, 1.);
                            acc_r.setValue(ijk, 1. );
                            acc_g.setValue(ijk, 1. );
                            acc_b.setValue(ijk, 1. );
                        }

*/
                        // acc_r.setValue(ijk, .0);
                        // acc_g.setValue(ijk, .0);
                        // acc_b.setValue(ijk, .0);
/*
                        // if ( interp_param_to_pixel.maxCoeff() > .9)
                        // {
                        //    acc_smoke_color.setValue(ijk, openvdb::Vec3s(float( interp_param_to_pixel(0) ), float( interp_param_to_pixel(1) ), float( interp_param_to_pixel(2) )) );
                            
                            if ( interp_param_to_pixel(1) < .1 && interp_param_to_pixel(2) < .1 )
                            {
                                acc_r.setValue(ijk, .9);
                                acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                            }
                            else if ( interp_param_to_pixel(0) < .1 && interp_param_to_pixel(2) < .1 )
                            {
                                acc_g.setValue(ijk, .9);
                                acc_strength.setValue(ijk, float( interp_param_to_pixel(1)) );

                            }
                            else if ( interp_param_to_pixel(0) < .1 && interp_param_to_pixel(1) < .1 )
                            {
                                acc_b.setValue(ijk, .9);
                                acc_strength.setValue(ijk, float( interp_param_to_pixel(2)) );
                            }
                        // }
*/
                    }

                }
            }
        }



    }


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