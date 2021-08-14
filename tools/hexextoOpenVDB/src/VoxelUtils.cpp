#include "VoxelUtils.h"
#include "TetMeshConnectivity.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

#include <openvdb/openvdb.h>


bool pointInsideT(const Eigen::Vector3d& A, 
                  const Eigen::Vector3d& B, 
                  const Eigen::Vector3d& C, 
                  const Eigen::Vector3d& D, 
                  const Eigen::Vector3d& samplePoint, Eigen::VectorXd& textureCoordinate)
{
    // Using this version: https://stackoverflow.com/a/51733522
    Eigen::MatrixXd Maff = Eigen::Matrix4d::Zero();
    Maff.block<1,3>(0,0) = B - A;
    Maff.block<1,3>(1,0) = C - A;
    Maff.block<1,3>(2,0) = D - A;
    Maff.block<1,3>(3,0) =     A;


// Double check how the other code does this.
    // Maff.block<1,3>(1,0) = B - A;
    // Maff.block<1,3>(2,0) = C - A;
    // Maff.block<1,3>(3,0) = D - A;
    // Maff.block<1,3>(0,0) =     A;

    Maff(3,3) = 1;

  //  std::cout << Maff << std::endl << std::endl;

    Eigen::MatrixXd M1 = Maff.transpose().inverse();

    Eigen::Vector4d p1 = Eigen::Vector4d::Constant(1.);
    p1.head(3) = samplePoint;

    Eigen::Vector4d paramPoint = M1 * p1;

    const Eigen::VectorXd& tC = paramPoint.head(3);

    bool isIn = (tC.maxCoeff() < 1.000 && tC.minCoeff() > -0.000);

    Eigen::Vector3d normal = (B-A).cross(C-A);

    isIn = isIn && ( tC.sum() < 1.0000); // tet has l1 constraint.  

    textureCoordinate = Eigen::VectorXd(4);
    textureCoordinate = paramPoint; 
    // std::cout << textureCoordinate.transpose() << std::endl;
    return isIn;

}

// weird thing that's kinda like a scaled jacobian.  Mostly for testing.  
Eigen::Matrix3d JacobianMetric(SceneInfo sc, int t_idx)
{
        // compute something like the scaled jacobian of the true map for testing  

        Eigen::Vector3d A_world = sc.V_unitcell.row( sc.T(t_idx, 0) );
        Eigen::Vector3d B_world = sc.V_unitcell.row( sc.T(t_idx, 1) );
        Eigen::Vector3d C_world = sc.V_unitcell.row( sc.T(t_idx, 2) );
        Eigen::Vector3d D_world = sc.V_unitcell.row( sc.T(t_idx, 3) );

        Eigen::Vector3d A_param = sc.param_unitcell.row( 4*t_idx + 0  );
        Eigen::Vector3d B_param = sc.param_unitcell.row( 4*t_idx + 1 );
        Eigen::Vector3d C_param = sc.param_unitcell.row( 4*t_idx + 2 );
        Eigen::Vector3d D_param = sc.param_unitcell.row( 4*t_idx + 3 );

        Eigen::Matrix3d toWorld;
        toWorld << B_world - A_world, C_world - A_world, D_world - A_world;

        Eigen::Matrix3d toParam;
        toParam << B_param - A_param, C_param - A_param, D_param - A_param;


        Eigen::Matrix3d jacobian = toParam * toWorld.inverse();

        return jacobian;
}


void stampScalarView(SceneInfo sc)
{
    double BIG_NUM = 100000000000000000.0;
    double line_w = .05;

    sc.param_pixel = sc.param_pixel * 1.; // rescale factor / max # of cells in a row.
    Eigen::MatrixXd param_gridcell = sc.param_unitcell * sc.cells;

    double shell_thick = .05;
    double border_w = .35;
    double facet_w = .01;
    double halo_w = .01;
    double cosmic_background = 0.0000;
    double scale_fac = .05;

    int ntets = sc.nTets();
    // int nverts = sc.nVerts();

    openvdb::FloatGrid::Accessor acc_r = sc.grid.r->getAccessor();
    openvdb::FloatGrid::Accessor acc_g = sc.grid.g->getAccessor();
    openvdb::FloatGrid::Accessor acc_b = sc.grid.b->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_r = sc.grid.smoke_r->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_g = sc.grid.smoke_g->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_b = sc.grid.smoke_b->getAccessor();
    openvdb::FloatGrid::Accessor acc_strength = sc.grid.strength->getAccessor(); // for now blackbody and emission are the same
    openvdb::FloatGrid::Accessor acc_smoke_density = sc.grid.smoke_density->getAccessor();
    openvdb::Vec3SGrid::Accessor acc_smoke_color = sc.grid.smoke_color->getAccessor();



    if( sc.V_curr == embedding::PARAM_SPACE )
    {
        std::cout << "plot param space" << std::endl;
    }
    else if( sc.V_curr == embedding::WORLD_SPACE )
    {
        std::cout << "plot world space" << std::endl;
    }




    for (int t = 0; t < ntets; t++)
    {
        // compute tet bbox 
        Eigen::Vector3d curt_min(BIG_NUM,BIG_NUM,BIG_NUM);
        Eigen::Vector3d curt_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);
        for (int v_idx = 0; v_idx < 4; v_idx++)
        {
            int rowId = sc.T(t, v_idx);
            // Eigen::Vector3d cur_point = sc.V_pixel.row( 4*t + v_idx );

            Eigen::Vector3d cur_point = sc.get_V_pos(t, v_idx);
            for ( int k = 0; k < 3; k++ )
            {
                if( cur_point(k) < curt_min(k) )
                    curt_min(k) = cur_point(k);
                if( cur_point(k) > curt_max(k) )
                    curt_max(k) = cur_point(k);
            }
        }


        // Eigen::Vector3d A = sc.V_pixel.row( sc.T(t, 0) );
        // Eigen::Vector3d B = sc.V_pixel.row( sc.T(t, 1) );
        // Eigen::Vector3d C = sc.V_pixel.row( sc.T(t, 2) );
        // Eigen::Vector3d D = sc.V_pixel.row( sc.T(t, 3) );


// // return vertex positions in pixel space.  
        Eigen::Vector3d A = sc.get_V_pos(t, 0 );
        Eigen::Vector3d B = sc.get_V_pos(t, 1 );
        Eigen::Vector3d C = sc.get_V_pos(t, 2 );
        Eigen::Vector3d D = sc.get_V_pos(t, 3 );

        Eigen::Matrix3d jacobian = JacobianMetric(sc, t);

        double jac0 = jacobian.col(0).norm();
        double jac1 = jacobian.col(1).norm();
        double jac2 = jacobian.col(2).norm();


        // double jac_min = std::min( std::min(jac0, jac1), jac2 ); 


        double jac_max = std::max( std::max(jac0, jac1), jac2 ); 

        jacobian = jacobian * ( 1. / jac_max );
        jac0 = jac0 / jac_max;
        jac1 = jac1 / jac_max;
        jac2 = jac2 / jac_max;


        double jac_det = jacobian.determinant(); 
        jac_det = std::sqrt( (1. - jac_det) * (1. - jac_det) );



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

                        acc_r.setValue(ijk, jac0 );
                        acc_g.setValue(ijk, jac1 );
                        acc_b.setValue(ijk, jac2 );
                        acc_strength.setValue(ijk, jac_det * jac_det * scale_fac);
         
                        acc_smoke_density.setValue(ijk, jac_det * jac_det * .001);
                        acc_smoke_r.setValue(ijk, jac0 );
                        acc_smoke_g.setValue(ijk, jac1 );
                        acc_smoke_b.setValue(ijk, jac2 );
//                         Eigen::Vector3d da = param_gridcell.row(4 * t + 1) - param_gridcell.row(4 * t);
//                         Eigen::Vector3d db = param_gridcell.row(4 * t + 2) - param_gridcell.row(4 * t);
//                         Eigen::Vector3d dc = param_gridcell.row(4 * t + 3) - param_gridcell.row(4 * t);
//                         Eigen::Vector3d orig = param_gridcell.row(4 * t);

//                         double s0 = param_val_tet(0);
//                         double s1 = param_val_tet(1);
//                         double s2 = param_val_tet(2);

//                         Eigen::Vector3d interp_param_to_pixel = s0 * da + s1 * db + s2 * dc + orig;
//                         interp_param_to_pixel *= 1. / 6.28;


//                         // This plots the param param per tet in the underlying mesh.  
//                         // Per tet normal vectors.  
//                         acc_smoke_color.setValue(ijk, openvdb::Vec3s(float( param_val_tet(0) ), 
//                                                                      float( param_val_tet(1) ), 
//                                                                      float( param_val_tet(2) )) );

                   
//                       //   std::cout << "interp_param_to_pixel" << interp_param_to_pixel.transpose() << std::endl<< std::endl;

//                         // Find texture coordinate of pixel rounded to a unit cube.  
//                         double unit_u = interp_param_to_pixel(0) - floor(interp_param_to_pixel(0));
//                         double unit_v = interp_param_to_pixel(1) - floor(interp_param_to_pixel(1));
//                         double unit_w = interp_param_to_pixel(2) - floor(interp_param_to_pixel(2));

// // Project vector onto closest line.  

//                         // multiply by jacobian. 

//                         double u_dist = std::min(unit_u, 1. - unit_u );
//                         double v_dist = std::min(unit_v, 1. - unit_v );
//                         double w_dist = std::min(unit_w, 1. - unit_w );






// ///////////////////////////////////////////////////////
//                         ///////////////////////  Plot GRID
//                         ///////////////////////////////

//                         // Stamp the colored grid pattern.
//                         if (sc.stamp_grid)
//                         {
//                             // acc_smoke_r.setValue(ijk, 0.9 );
//                             // acc_smoke_g.setValue(ijk, 0.9  );
//                             // acc_smoke_b.setValue(ijk, 0.9  );
//                             // acc_smoke_density.setValue(ijk, 0.000 );

//                             if ( v_dist < line_w && w_dist < line_w )
//                             {
//                                 acc_r.setValue(ijk, .6 );
//                                 acc_g.setValue(ijk, .0 );
//                                 acc_b.setValue(ijk, 0. );
//                                 acc_strength.setValue(ijk, unit_u);
         
//                                 acc_smoke_density.setValue(ijk, 1.);
//                                 acc_smoke_r.setValue(ijk, 1. );
//                                 acc_smoke_g.setValue(ijk, .9 );
//                                 acc_smoke_b.setValue(ijk, .9 );
//                                 // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
//                             }
//                             else if ( u_dist < line_w && w_dist < line_w )
//                             {
//                                 acc_r.setValue(ijk, .0 );
//                                 acc_g.setValue(ijk, .6 );
//                                 acc_b.setValue(ijk, 0. );

//                                 acc_strength.setValue(ijk, unit_v);

//                                 acc_smoke_density.setValue(ijk, 1.);
//                                 acc_smoke_r.setValue(ijk, .9 );
//                                 acc_smoke_g.setValue(ijk, 1. );
//                                 acc_smoke_b.setValue(ijk, .9 );

//                             }
//                             else if ( u_dist < line_w && v_dist < line_w )
//                             {
//                                 acc_r.setValue(ijk, .0 );
//                                 acc_g.setValue(ijk, .0 );
//                                 acc_b.setValue(ijk, .6 );

//                                 acc_strength.setValue(ijk, unit_w);

//                                 acc_smoke_density.setValue(ijk, 1.);
//                                 acc_smoke_r.setValue(ijk, .9 );
//                                 acc_smoke_g.setValue(ijk, .9 );
//                                 acc_smoke_b.setValue(ijk, 1.0 );
//                             }


//                         }

                     }
                }
            }
        }
    }
}





void stampParamView(SceneInfo sc)
{
    double BIG_NUM = 100000000000000000.0;
    double line_w = .05;

    sc.param_pixel = sc.param_pixel * 1.; // rescale factor / max # of cells in a row.
    Eigen::MatrixXd param_gridcell = sc.param_unitcell * sc.cells;


    double shell_thick = .05;
    double space_thick = shell_thick / 2.;

    double border_w = .4;
    double facet_w = .01;
    double halo_w = .01;
    double cosmic_background = 0.0000;

    int ntets = sc.nTets();
    // int nverts = sc.nVerts();

    openvdb::FloatGrid::Accessor acc_r = sc.grid.r->getAccessor();
    openvdb::FloatGrid::Accessor acc_g = sc.grid.g->getAccessor();
    openvdb::FloatGrid::Accessor acc_b = sc.grid.b->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_r = sc.grid.smoke_r->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_g = sc.grid.smoke_g->getAccessor();
    openvdb::FloatGrid::Accessor acc_smoke_b = sc.grid.smoke_b->getAccessor();
    openvdb::FloatGrid::Accessor acc_strength = sc.grid.strength->getAccessor(); // for now blackbody and emission are the same
    openvdb::FloatGrid::Accessor acc_smoke_density = sc.grid.smoke_density->getAccessor();
    openvdb::Vec3SGrid::Accessor acc_smoke_color = sc.grid.smoke_color->getAccessor();



    if( sc.V_curr == embedding::PARAM_SPACE )
    {
        std::cout << "plot param space" << std::endl;
    }
    else if( sc.V_curr == embedding::WORLD_SPACE )
    {
        std::cout << "plot world space" << std::endl;
    }




    for (int t = 0; t < ntets; t++)
    {
        // compute tet bbox 
        Eigen::Vector3d curt_min(BIG_NUM,BIG_NUM,BIG_NUM);
        Eigen::Vector3d curt_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);
        for (int v_idx = 0; v_idx < 4; v_idx++)
        {
            int rowId = sc.T(t, v_idx);
            // Eigen::Vector3d cur_point = sc.V_pixel.row( 4*t + v_idx );

            Eigen::Vector3d cur_point = sc.get_V_pos(t, v_idx);
            for ( int k = 0; k < 3; k++ )
            {
                if( cur_point(k) < curt_min(k) )
                    curt_min(k) = cur_point(k);
                if( cur_point(k) > curt_max(k) )
                    curt_max(k) = cur_point(k);
            }
        }


        // Eigen::Vector3d A = sc.V_pixel.row( sc.T(t, 0) );
        // Eigen::Vector3d B = sc.V_pixel.row( sc.T(t, 1) );
        // Eigen::Vector3d C = sc.V_pixel.row( sc.T(t, 2) );
        // Eigen::Vector3d D = sc.V_pixel.row( sc.T(t, 3) );


// // return vertex positions in pixel space.  
        Eigen::Vector3d A = sc.get_V_pos(t, 0 );
        Eigen::Vector3d B = sc.get_V_pos(t, 1 );
        Eigen::Vector3d C = sc.get_V_pos(t, 2 );
        Eigen::Vector3d D = sc.get_V_pos(t, 3 );




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
               
                        Eigen::Vector3d da = param_gridcell.row(4 * t + 1) - param_gridcell.row(4 * t);
                        Eigen::Vector3d db = param_gridcell.row(4 * t + 2) - param_gridcell.row(4 * t);
                        Eigen::Vector3d dc = param_gridcell.row(4 * t + 3) - param_gridcell.row(4 * t);
                        Eigen::Vector3d orig = param_gridcell.row(4 * t);

                        double s0 = param_val_tet(0);
                        double s1 = param_val_tet(1);
                        double s2 = param_val_tet(2);

                        Eigen::Vector3d interp_param_to_pixel = s0 * da + s1 * db + s2 * dc + orig;
                        interp_param_to_pixel *= 1. / 6.28;


                        // This plots the param param per tet in the underlying mesh.  
                        // Per tet normal vectors.  
                        acc_smoke_color.setValue(ijk, openvdb::Vec3s(float( param_val_tet(0) ), 
                                                                     float( param_val_tet(1) ), 
                                                                     float( param_val_tet(2) )) );

                   
                      //   std::cout << "interp_param_to_pixel" << interp_param_to_pixel.transpose() << std::endl<< std::endl;

                        // Find texture coordinate of pixel rounded to a unit cube.  
                        double unit_u = interp_param_to_pixel(0) - floor(interp_param_to_pixel(0));
                        double unit_v = interp_param_to_pixel(1) - floor(interp_param_to_pixel(1));
                        double unit_w = interp_param_to_pixel(2) - floor(interp_param_to_pixel(2));

// Project vector onto closest line.  

                        // multiply by jacobian. 

                        double u_dist = std::min(unit_u, 1. - unit_u );
                        double v_dist = std::min(unit_v, 1. - unit_v );
                        double w_dist = std::min(unit_w, 1. - unit_w );






///////////////////////////////////////////////////////
                        ///////////////////////  Plot GRID
                        ///////////////////////////////

                        // Stamp the colored grid pattern.
                        if (sc.stamp_grid)
                        {
                            // acc_smoke_r.setValue(ijk, 0.9 );
                            // acc_smoke_g.setValue(ijk, 0.9  );
                            // acc_smoke_b.setValue(ijk, 0.9  );
                            // acc_smoke_density.setValue(ijk, 0.000 );

                            if ( v_dist < line_w && w_dist < line_w )
                            {
                                acc_r.setValue(ijk, .6 );
                                acc_g.setValue(ijk, .0 );
                                acc_b.setValue(ijk, 0. );
                                acc_strength.setValue(ijk, unit_u);
         
                                acc_smoke_density.setValue(ijk, 1.);
                                acc_smoke_r.setValue(ijk, 1. );
                                acc_smoke_g.setValue(ijk, .9 );
                                acc_smoke_b.setValue(ijk, .9 );
                                // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                            }
                            else if ( u_dist < line_w && w_dist < line_w )
                            {
                                acc_r.setValue(ijk, .0 );
                                acc_g.setValue(ijk, .6 );
                                acc_b.setValue(ijk, 0. );

                                acc_strength.setValue(ijk, unit_v);

                                acc_smoke_density.setValue(ijk, 1.);
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, 1. );
                                acc_smoke_b.setValue(ijk, .9 );

                            }
                            else if ( u_dist < line_w && v_dist < line_w )
                            {
                                acc_r.setValue(ijk, .0 );
                                acc_g.setValue(ijk, .0 );
                                acc_b.setValue(ijk, .6 );

                                acc_strength.setValue(ijk, unit_w);

                                acc_smoke_density.setValue(ijk, 1.);
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, .9 );
                                acc_smoke_b.setValue(ijk, 1.0 );
                            }


                        }


///////////////////////////////////////////////////////
                        ///////////////////////  Plot Exploded Hexes
                        ///////////////////////////////

                        if (sc.stamp_centers)
                        {


                            if (u_dist > border_w && v_dist > border_w && w_dist > border_w)
                            {
                                acc_smoke_r.setValue(ijk, 0.85 );
                                acc_smoke_g.setValue(ijk, 0.85  );
                                acc_smoke_b.setValue(ijk, 0.85  );
                                acc_smoke_density.setValue(ijk, .6 );

                                acc_r.setValue(ijk, .55 );
                                acc_g.setValue(ijk, .55 );
                                acc_b.setValue(ijk, .55 );
                                acc_strength.setValue(ijk, .1);
                            }
                            else if ( u_dist > border_w - space_thick && 
                                v_dist > border_w && 
                                w_dist > border_w)
                            {
                                acc_r.setValue(ijk, .0 );
                                acc_g.setValue(ijk, .7 );
                                acc_b.setValue(ijk, .7 );
                                acc_strength.setValue(ijk, .0);

                                acc_smoke_r.setValue(ijk, 0.95 );
                                acc_smoke_g.setValue(ijk, 1.0  );
                                acc_smoke_b.setValue(ijk, 1.  );
                                acc_smoke_density.setValue(ijk, .0 );
                            }
                            else if ( u_dist > border_w && 
                                      v_dist > border_w - space_thick && 
                                      w_dist > border_w)
                            {
                                acc_r.setValue(ijk, .7);
                                acc_g.setValue(ijk, .0 );
                                acc_b.setValue(ijk, .7 );
                                acc_strength.setValue(ijk, 0);

                                acc_smoke_r.setValue(ijk, 1.0 );
                                acc_smoke_g.setValue(ijk, 0.95  );
                                acc_smoke_b.setValue(ijk, 1.  );
                                acc_smoke_density.setValue(ijk, .00);
                            }
                            else if ( u_dist > border_w && 
                                      v_dist > border_w && 
                                      w_dist > border_w - space_thick)
                            {
                                acc_r.setValue(ijk, .7 );
                                acc_g.setValue(ijk, .7 );
                                acc_b.setValue(ijk, .0 );
                                acc_strength.setValue(ijk, .0);

                                acc_smoke_r.setValue(ijk, 1.0 );
                                acc_smoke_g.setValue(ijk, 1.0  );
                                acc_smoke_b.setValue(ijk, .95  );
                                acc_smoke_density.setValue(ijk, .0 );
                            }
                            else if ( u_dist > border_w - shell_thick && 
                                      v_dist > border_w && 
                                      w_dist > border_w)
                            {
                                acc_r.setValue(ijk, .0 );
                                acc_g.setValue(ijk, .7 );
                                acc_b.setValue(ijk, .7 );
                                acc_strength.setValue(ijk, .6);

                                acc_smoke_r.setValue(ijk, 0.5 );
                                acc_smoke_g.setValue(ijk, .9  );
                                acc_smoke_b.setValue(ijk, .9  );
                                acc_smoke_density.setValue(ijk, .1 );
                            }
                            else if ( u_dist > border_w && 
                                      v_dist > border_w - shell_thick && 
                                      w_dist > border_w)
                            {
                                acc_r.setValue(ijk, .7);
                                acc_g.setValue(ijk, .0 );
                                acc_b.setValue(ijk, .7 );
                                acc_strength.setValue(ijk, .6);

                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, 0.5 );
                                acc_smoke_b.setValue(ijk, .9  );
                                acc_smoke_density.setValue(ijk, .1 );
                            }
                            else if ( u_dist > border_w && 
                                      v_dist > border_w && 
                                      w_dist > border_w - shell_thick)
                            {
                                acc_r.setValue(ijk, .7 );
                                acc_g.setValue(ijk, .7 );
                                acc_b.setValue(ijk, .0 );
                                acc_strength.setValue(ijk, .6);

                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, .9  );
                                acc_smoke_b.setValue(ijk, .5  );
                                acc_smoke_density.setValue(ijk, .1 );
                            }






                        }


                    }

                }
            }
        }
    }


}

