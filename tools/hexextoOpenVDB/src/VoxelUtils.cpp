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

    bool isIn = (tC.maxCoeff() < 1.000 && tC.minCoeff() > -0.001);

    // Eigen::Vector3d normal = (B-A).cross(C-A);

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
        // Eigen::Vector3d A_param = sc.param_unitcell.row( 4*t_idx + 3  );
        // Eigen::Vector3d B_param = sc.param_unitcell.row( 4*t_idx + 2 );
        // Eigen::Vector3d C_param = sc.param_unitcell.row( 4*t_idx + 1 );
        // Eigen::Vector3d D_param = sc.param_unitcell.row( 4*t_idx + 0 );


        Eigen::Matrix3d toWorld;
        toWorld << B_world - A_world, C_world - A_world, D_world - A_world;
        // toWorld = toWorld.transpose();


        // std::cout << "world " << toWorld << std::endl;

        // std::cout << "world row1 " << B_world - A_world << std::endl;
        // std::cout << "world row1 " << C_world - A_world << std::endl;
        // std::cout << "world row1 " << D_world - A_world << std::endl;


        Eigen::Matrix3d toParam;
        toParam << B_param - A_param, C_param - A_param, D_param - A_param;
        // toParam = toParam.transpose();


        Eigen::Matrix3d jacobian =  toParam * toWorld.inverse();

        return jacobian;
}

/*

extractSingularCurveNetwork(const Eigen::MatrixXd& V, 
    const CubeCover::TetMeshConnectivity& mesh, 
    const CubeCover::FrameField& field, 
    Eigen::MatrixXd& Pgreen, Eigen::MatrixXi& Egreen,
    Eigen::MatrixXd& Pblue, Eigen::MatrixXi& Eblue,
    Eigen::MatrixXd& Pblack, Eigen::MatrixXi& Eblack
)

*/


//////////////////////\
/////\  Stamp Surface
////////////////////

void stampSurfaceView(SceneInfo sc, CubeCover::TetMeshConnectivity mesh, CubeCover::FrameField* field)
{
    double BIG_NUM = 100000000000000000.0;

    sc.param_pixel = sc.param_pixel * 1.; // rescale factor / max # of cells in a row.
    Eigen::MatrixXd param_gridcell = sc.param_unitcell * sc.cells;


    double scale_fac = .05;
    double smoke_scale_fac = .0005;

    int thickness = 10;

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

    // int tris = bdryF.rows();




    for (int t = 0; t < ntets; t++)
    {
        for (int idx = 0; idx < 4; idx++)
        {
            int curFace = mesh.tetFace(t, idx);
            if (!field->faceAssignment(curFace).isIdentity())
            {
                int tidx = 0;
                if (mesh.faceTet(curFace, 1) == t)
                {
                    tidx = 1;
                }

                int v0 = mesh.faceTetVertexIndex(curFace, tidx, 0);
                int v1 = mesh.faceTetVertexIndex(curFace, tidx, 1);
                int v2 = mesh.faceTetVertexIndex(curFace, tidx, 2);

                Eigen::Vector3d A = sc.get_V_pos(t, v0 );
                Eigen::Vector3d B = sc.get_V_pos(t, v1 );
                Eigen::Vector3d C = sc.get_V_pos(t, v2 );

                if ( mesh.faceTet(curFace, 0) == -1 )
                {
                    Eigen::Vector3d tmp = B;
                    B = A; 
                    A = tmp;
                }

                int xsamp = floor( (B - A).norm() );
                int ysamp = floor( (C - A).norm() );

                Eigen::Vector3d face_base_1 = (B - A) / (B - A).norm();
                Eigen::Vector3d face_base_2 = (C - A) / (C - A).norm();
                Eigen::Vector3d face_n = face_base_1.cross(face_base_2);


                openvdb::Coord ijk;
                for (int i = 0; i < xsamp; i++)
                {
                    for (int j = 0; j < ysamp - floor( ysamp * double(i) / double(xsamp) ); j++)
                    {
                        for (int k = 0; k < thickness; k++)
                        {

                            Eigen::Vector3d pos = A + i * face_base_1
                                                    + j * face_base_2
                                                    + k * face_n;


                            ijk[0] = floor(pos[0]);
                            ijk[1] = floor(pos[1]);
                            ijk[2] = floor(pos[2]);


                            double param_tol = .01;
                            if ( tidx == 0)
                            {
                                acc_r.setValue(ijk, .0 );
                                acc_g.setValue(ijk, .5 );
                                acc_b.setValue(ijk, .9 );
                                acc_strength.setValue(ijk, 0.2);
             
                                acc_smoke_density.setValue(ijk, .1);
                                acc_smoke_r.setValue(ijk, .0 );
                                acc_smoke_g.setValue(ijk, 0.2 );
                                acc_smoke_b.setValue(ijk, 0.5 );
                            }
                            else if ( tidx == 1 )
                            {
                                acc_r.setValue(ijk, .9 );
                                acc_g.setValue(ijk, .5 );
                                acc_b.setValue(ijk, .0 );
                                acc_strength.setValue(ijk, 0.2);
             
                                acc_smoke_density.setValue(ijk, .1);
                                acc_smoke_r.setValue(ijk, .5 );
                                acc_smoke_g.setValue(ijk, 0.2 );
                                acc_smoke_b.setValue(ijk, 0.0 );
                            }

     



                        }

                    }
                }

            }


            if ( mesh.isBoundaryFace(curFace) )
            {
                // first, just plot the face.  

                // Eigen::Vector3d A = sc.get_V_pos(t, 0 );
                // Eigen::Vector3d B = sc.get_V_pos(t, 1 );
                // Eigen::Vector3d C = sc.get_V_pos(t, 2 );

                int tidx = 0;
                if ( mesh.faceTet(curFace, 0) == -1 )
                    tidx = 1;

                int v0 = mesh.faceTetVertexIndex(curFace, tidx, 0);
                int v1 = mesh.faceTetVertexIndex(curFace, tidx, 1);
                int v2 = mesh.faceTetVertexIndex(curFace, tidx, 2);

                Eigen::Vector3d A = sc.get_V_pos(t, v0 );
                Eigen::Vector3d B = sc.get_V_pos(t, v1 );
                Eigen::Vector3d C = sc.get_V_pos(t, v2 );

                if ( mesh.faceTet(curFace, 0) == -1 )
                {
                    Eigen::Vector3d tmp = B;
                    B = A; 
                    A = tmp;
                }

                int xsamp = floor( (B - A).norm() );
                int ysamp = floor( (C - A).norm() );

                Eigen::Vector3d face_base_1 = (B - A) / (B - A).norm();
                Eigen::Vector3d face_base_2 = (C - A) / (C - A).norm();
                Eigen::Vector3d face_n = face_base_1.cross(face_base_2);

                Eigen::Vector3d pa = param_gridcell.row(4 * t + v0);
                Eigen::Vector3d pb = param_gridcell.row(4 * t + v1);
                Eigen::Vector3d pc = param_gridcell.row(4 * t + v2);

                Eigen::Vector3d du = pb - pa; 
                Eigen::Vector3d dv = pc - pa;
                du = du.cwiseAbs();
                dv = dv.cwiseAbs();


                openvdb::Coord ijk;
                for (int i = 0; i < xsamp; i++)
                {
                    for (int j = 0; j < ysamp - floor( ysamp * double(i) / double(xsamp) ); j++)
                    {
                        for (int k = 0; k < thickness; k++)
                        {

                            Eigen::Vector3d pos = A + i * face_base_1
                                                    + j * face_base_2
                                                    + k * face_n;


                            ijk[0] = floor(pos[0]);
                            ijk[1] = floor(pos[1]);
                            ijk[2] = floor(pos[2]);


                            double param_tol = .01;
                            if ( du(0) < param_tol && dv(0) < param_tol)
                            {
                                acc_r.setValue(ijk, .2 );
                                acc_g.setValue(ijk, .9 );
                                acc_b.setValue(ijk, .9 );
                                acc_strength.setValue(ijk, 0.03);
             
                                acc_smoke_density.setValue(ijk, .01);
                                acc_smoke_r.setValue(ijk, .0 );
                                acc_smoke_g.setValue(ijk, 0.9 );
                                acc_smoke_b.setValue(ijk, 0.9 );
                            }
                            else if ( du(1) < param_tol && dv(1) < param_tol)
                            {
                                acc_r.setValue(ijk, .9 );
                                acc_g.setValue(ijk, .2 );
                                acc_b.setValue(ijk, .9 );
                                acc_strength.setValue(ijk, 0.03);
             
                                acc_smoke_density.setValue(ijk, .01);
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, 0.0 );
                                acc_smoke_b.setValue(ijk, 0.9 );
                            }

                            else if ( du(2) < param_tol && dv(2) < param_tol)
                            {
                                acc_r.setValue(ijk, .9 );
                                acc_g.setValue(ijk, .9 );
                                acc_b.setValue(ijk, .2 );
                                acc_strength.setValue(ijk, 0.03);
             
                                acc_smoke_density.setValue(ijk, .01);
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, 0.9 );
                                acc_smoke_b.setValue(ijk, 0. );
                            }



                        }

                    }
                }


            }
        }


    }

}




void stampScalarView(SceneInfo sc)
{
    double BIG_NUM = 100000000000000000.0;

    sc.param_pixel = sc.param_pixel * 1.; // rescale factor / max # of cells in a row.
    Eigen::MatrixXd param_gridcell = sc.param_unitcell * sc.cells;


    double scale_fac = .05;
    double smoke_scale_fac = .025;

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
         
                        acc_smoke_density.setValue(ijk, jac_det * jac_det * smoke_scale_fac);
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
    

    sc.param_pixel = sc.param_pixel * 1.; // rescale factor / max # of cells in a row.
    Eigen::MatrixXd param_gridcell = sc.param_unitcell * sc.cells;


    double shell_thick = .05;
    double space_thick = shell_thick / 2.;
    double line_w = .05;
    double world_line_w = .05;

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



    std::cout << " number of tets in mesh " << ntets << std::endl;
     for (int t = 0; t < ntets; t++)
     {
   // int t = 17;
        bool tetIsActive = false;

        // compute tet bbox 
        Eigen::Vector3d curt_min(BIG_NUM,BIG_NUM,BIG_NUM);
        Eigen::Vector3d curt_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);
        for (int v_idx = 0; v_idx < 4; v_idx++)
        {
          //  int rowId = sc.T(t, v_idx);
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

        // std::cout << curt_min << "cur t_min" << std::endl;
        // std::cout << curt_max << "cur t_max" << std::endl;
        


        // Eigen::Vector3d A = sc.V_pixel.row( sc.T(t, 0) );
        // Eigen::Vector3d B = sc.V_pixel.row( sc.T(t, 1) );
        // Eigen::Vector3d C = sc.V_pixel.row( sc.T(t, 2) );
        // Eigen::Vector3d D = sc.V_pixel.row( sc.T(t, 3) );


// // return vertex positions in pixel space.  
        Eigen::Vector3d A = sc.get_V_pos(t, 0 );
        Eigen::Vector3d B = sc.get_V_pos(t, 1 );
        Eigen::Vector3d C = sc.get_V_pos(t, 2 );
        Eigen::Vector3d D = sc.get_V_pos(t, 3 );

        double det = (B-A).cross(C-A).dot(D-A);

        if (det < .0000005)
        { 
            Eigen::Vector3d tmp = B;
            B = C; 
            C = tmp;

            // std::cout << " Small tet determinant " << det << std::endl;
            // curt_min = Eigen::Vector3d::Zero();
            // curt_max = Eigen::Vector3d::Zero();
        }

                    Eigen::Vector3d tmp = B;
            B = C; 
            C = tmp;

  
        // if (det < .0000005)
        // {      
        //     // std::cout << " Small tet determinant " << det << std::endl;
        //     curt_min = Eigen::Vector3d::Zero();
        //     curt_max = Eigen::Vector3d::Zero();
        // }


        Eigen::Matrix3d jacobian = JacobianMetric(sc, t);  // J * Param  = World  

        // this is used to rescale the sampling.  

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
                    // if ( pIsIn )
                    // {
               
                        Eigen::Vector3d da = param_gridcell.row(4 * t + 1) - param_gridcell.row(4 * t);
                        Eigen::Vector3d db = param_gridcell.row(4 * t + 2) - param_gridcell.row(4 * t);
                        Eigen::Vector3d dc = param_gridcell.row(4 * t + 3) - param_gridcell.row(4 * t);
                        Eigen::Vector3d orig = param_gridcell.row(4 * t);

                        double s0 = param_val_tet(0);
                        double s1 = param_val_tet(1);
                        double s2 = param_val_tet(2);

                        Eigen::Vector3d interp_param_to_pixel = s0 * da + s1 * db + s2 * dc + orig;
                        interp_param_to_pixel *= 1. / 6.28;


                        // This plots the param per tet in the underlying mesh.  
                        // Per tet normal vectors.  
                        // acc_smoke_color.setValue(ijk, openvdb::Vec3s(float( param_val_tet(0) ), 
                        //                                              float( param_val_tet(1) ), 
                        //                                              float( param_val_tet(2) )) );

                   
                      //   std::cout << "interp_param_to_pixel" << interp_param_to_pixel.transpose() << std::endl<< std::endl;

                        // Find texture coordinate of pixel rounded to a unit cube.  
                        double unit_u = interp_param_to_pixel(0) - floor(interp_param_to_pixel(0));
                        double unit_v = interp_param_to_pixel(1) - floor(interp_param_to_pixel(1));
                        double unit_w = interp_param_to_pixel(2) - floor(interp_param_to_pixel(2));

// Project vector onto closest line.  

                        // multiply by jacobian. 
                        // if (  ) 
                        double u_dist = std::min(unit_u, 1. - unit_u );
                        double v_dist = std::min(unit_v, 1. - unit_v );
                        double w_dist = std::min(unit_w, 1. - unit_w );

                        Eigen::Vector3d paramVal_shifted;
                        paramVal_shifted << u_dist, v_dist, w_dist;

                        Eigen::Vector3d param2world = jacobian * paramVal_shifted;

                        // if( sc.V_curr == embedding::WORLD_SPACE )
                        // {

                        //     u_dist = param2world[0];
                        //     v_dist = param2world[1];
                        //     w_dist = param2world[2];

                        // }



                        Eigen::Vector3d roundedWorld = paramVal_shifted;
                        if ( v_dist < line_w && w_dist < line_w )
                        {
                            roundedWorld[0] = 0.;
                        }
                        if ( u_dist < line_w && w_dist < line_w )
                        {
                            roundedWorld[1] = 0.;
                        }
                        if ( u_dist < line_w && v_dist < line_w )
                        {
                            roundedWorld[2] = 0.;
                        }
                        

                        // std::cout << paramVal_shifted << " param val shift " << std::endl;
                        // // std::cout << paramVal << " param val " << std::endl;
                         // std::cout << u_dist << " " << v_dist << " " << w_dist << " u, v, w dist " << std::endl;
                        // std::cout << param2world << " param2world" << std::endl;
                        // std::cout << roundedWorld << " roundedWorld " << std::endl;


                        if (sc.stamp_grid)
                        {

                                    // if ( pIsIn )
                                    // {
                                    //     acc_r.setValue(ijk, u_dist );
                                    //     acc_g.setValue(ijk, v_dist );
                                    //     acc_b.setValue(ijk, w_dist );
                                    // // acc_strength.setValue(ijk, unit_u);
                                    //     acc_strength.setValue(ijk, .05);

                                    //     acc_smoke_density.setValue(ijk, .1);
                                    //     acc_smoke_r.setValue(ijk, u_dist );
                                    //     acc_smoke_g.setValue(ijk, v_dist );
                                    //     acc_smoke_b.setValue(ijk, w_dist );



                                    // }
                             

                     //        if ( pIsIn ) // pIsIn && t % 17 == 0 )
                     //        {
                     //            acc_r.setValue(ijk, 1. );
                     //            acc_g.setValue(ijk, 0. );
                     //            acc_b.setValue(ijk, 0. );
                     //                // acc_strength.setValue(ijk, unit_u);
                     //            acc_strength.setValue(ijk, .2);

                     //            acc_smoke_density.setValue(ijk, .1);
                     //            acc_smoke_r.setValue(ijk, 1. );
                     //            acc_smoke_g.setValue(ijk, 1. );
                     //            acc_smoke_b.setValue(ijk, 1. );

                     // //           std::cout << ijk << " ijk " <<  "is in" << std::endl;



                     //        }
                     //        if ( acc_r.getValue(ijk) < .5 ) // t % 17 == 0 )
                     //        {
                     //            acc_r.setValue(ijk, 0. );
                     //            acc_g.setValue(ijk, 0. );
                     //            acc_b.setValue(ijk, .5 );
                     //            //     // acc_strength.setValue(ijk, unit_u);
                     //            acc_strength.setValue(ijk, .05);

                     //            acc_smoke_density.setValue(ijk, .1);
                     //            acc_smoke_r.setValue(ijk, .9 );
                     //            acc_smoke_g.setValue(ijk, .9 );
                     //            acc_smoke_b.setValue(ijk, .5 );

                     //    //        std::cout << ijk << " ijk " <<  "is out" << std::endl;
                     //        }

                            
                                if ( v_dist < line_w && w_dist < line_w )
                                {
                                    if ( pIsIn )
                                    {
                                        // tetIsActive = true;
                                        acc_r.setValue(ijk, .6 );
                                        acc_g.setValue(ijk, .0 );
                                        acc_b.setValue(ijk, 0. );

                                        // acc_r.setValue(ijk, .0 );
                                        // acc_g.setValue(ijk, .0 );
                                        // acc_b.setValue(ijk, 0.6 );
                                    // // acc_strength.setValue(ijk, unit_u);
                                    //     acc_strength.setValue(ijk, .5);

                                    //     acc_r.setValue(ijk, jac0 );
                                    //     acc_g.setValue(ijk, jac1 );
                                    //     acc_b.setValue(ijk, jac2 );
                                    // // acc_strength.setValue(ijk, unit_u);
                                        acc_strength.setValue(ijk, .05);

                                        acc_smoke_density.setValue(ijk, .1);
                                        acc_smoke_r.setValue(ijk, u_dist );
                                        acc_smoke_g.setValue(ijk, v_dist );
                                        acc_smoke_b.setValue(ijk, w_dist );


                                    }
                                    // else if ( acc_r.getValue(ijk) > 0.5 )
                                    // {
                                    //     Eigen::MatrixXd Maff = Eigen::Matrix4d::Zero();
                                    //     Maff.block<1,3>(0,0) = B - A;
                                    //     Maff.block<1,3>(1,0) = C - A;
                                    //     Maff.block<1,3>(2,0) = D - A;
                                    //     Maff.block<1,3>(3,0) =     A;


                                    //   //  acc_b.setValue(ijk, 0.6 );
                                    // // Double check how the other code does this.
                                    //     // Maff.block<1,3>(1,0) = B - A;
                                    //     // Maff.block<1,3>(2,0) = C - A;
                                    //     // Maff.block<1,3>(3,0) = D - A;
                                    //     // Maff.block<1,3>(0,0) =     A;

                                    //     Maff(3,3) = 1;

                                    //   //  std::cout << Maff << std::endl << std::endl;

                                    //     Eigen::MatrixXd M1 = Maff.transpose().inverse();

                                    //     Eigen::Vector4d p1 = Eigen::Vector4d::Constant(1.);
                                    //     p1.head(3) = p;

                                    //     Eigen::Vector4d paramPoint = M1 * p1;

                                    //     const Eigen::VectorXd& tC = paramPoint.head(3);

                                    //     bool isIn = (tC.maxCoeff() < 1.000 && tC.minCoeff() > -0.001);

                                    //     // Eigen::Vector3d normal = (B-A).cross(C-A);

                                    //     isIn = isIn && ( tC.sum() < 1.0000); // tet has l1 constraint.  
                                    //     std::cout << tC.maxCoeff() << " max " << tC.minCoeff() << " min " << tC.sum() << " sum " << std::endl;

                                    //     Eigen::VectorXd textureCoordinate = Eigen::VectorXd(4);
                                    //     textureCoordinate = paramPoint; 
                                    //     std::cout << " tex_coord " << textureCoordinate.transpose() << std::endl;
                                    // }
                                    else if ( acc_r.getValue(ijk) < .5 )
                                    {
                                        // acc_r.setValue(ijk, .0 );
                                        // acc_g.setValue(ijk, .0 );
                                        acc_b.setValue(ijk, 0.6 );
                                    // acc_strength.setValue(ijk, unit_u);
                                        acc_strength.setValue(ijk, .05);
                                    }

                                    if ( roundedWorld.norm() < world_line_w )
                                    {
                                        acc_smoke_density.setValue(ijk, 1.);
                                        acc_smoke_r.setValue(ijk, 1. );
                                        acc_smoke_g.setValue(ijk, .9 );
                                        acc_smoke_b.setValue(ijk, .9 );
                                    }
             
                                   
                                    // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                                }


                            
                            

       

        
                            // else if ( u_dist < line_w && w_dist < line_w )
                            // {
                            //     acc_r.setValue(ijk, .0 );
                            //     acc_g.setValue(ijk, .6 );
                            //     acc_b.setValue(ijk, 0. );

                            //     acc_strength.setValue(ijk, unit_v);

                            //     if ( roundedWorld.norm() < world_line_w )
                            //     {
                            //         acc_smoke_density.setValue(ijk, 1.);
                            //         acc_smoke_r.setValue(ijk, .9 );
                            //         acc_smoke_g.setValue(ijk, 1. );
                            //         acc_smoke_b.setValue(ijk, .9 );
                            //     }



                            // }
                            // else if ( u_dist < line_w && v_dist < line_w )
                            // {
                            //     acc_r.setValue(ijk, .0 );
                            //     acc_g.setValue(ijk, .0 );
                            //     acc_b.setValue(ijk, .6 );

                            //     acc_strength.setValue(ijk, unit_w);

                            //     if ( roundedWorld.norm() < world_line_w )
                            //     {
                            //         acc_smoke_density.setValue(ijk, 1.);
                            //         acc_smoke_r.setValue(ijk, .9 );
                            //         acc_smoke_g.setValue(ijk, .9 );
                            //         acc_smoke_b.setValue(ijk, 1.0 );
                            //     }


                            // }
                        }
               

                        // if (sc.stamp_grid)
                        // {
                        //    char color = 'r';
                        //    if ()
                        //    if( param2world.norm() < line_w ) 
                        //    {


                        //    }
                        // }



// test vector in the axis plane so it's 2D, evaluate jacobian to put into world coordiantes before filtering.





// ///////////////////////////////////////////////////////
//                         ///////////////////////  Plot GRID
//                         ///////////////////////////////

//                         // // Stamp the colored grid pattern.
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

//                                 param2world[0] = 0.;

//                                 if ( param2world.norm() < line_w )
//                                 {
//                                     acc_smoke_density.setValue(ijk, 1.);
//                                     acc_smoke_r.setValue(ijk, 1. );
//                                     acc_smoke_g.setValue(ijk, .9 );
//                                     acc_smoke_b.setValue(ijk, .9 );
//                                 }
         
                               
//                                 // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
//                             }
//                             else if ( u_dist < line_w && w_dist < line_w )
//                             {
//                                 acc_r.setValue(ijk, .0 );
//                                 acc_g.setValue(ijk, .6 );
//                                 acc_b.setValue(ijk, 0. );

//                                 acc_strength.setValue(ijk, unit_v);

//                                 param2world[1] = 0.;

//                                 if ( param2world.norm() < line_w )
//                                 {
//                                     acc_smoke_density.setValue(ijk, 1.);
//                                     acc_smoke_r.setValue(ijk, .9 );
//                                     acc_smoke_g.setValue(ijk, 1. );
//                                     acc_smoke_b.setValue(ijk, .9 );
//                                 }



//                             }
//                             else if ( u_dist < line_w && v_dist < line_w )
//                             {
//                                 acc_r.setValue(ijk, .0 );
//                                 acc_g.setValue(ijk, .0 );
//                                 acc_b.setValue(ijk, .6 );

//                                 acc_strength.setValue(ijk, unit_w);

//                                 param2world[2] = 0.;

//                                 if ( param2world.norm() < line_w )
//                                 {
//                                     acc_smoke_density.setValue(ijk, 1.);
//                                     acc_smoke_r.setValue(ijk, .9 );
//                                     acc_smoke_g.setValue(ijk, .9 );
//                                     acc_smoke_b.setValue(ijk, 1.0 );
//                                 }


//                             }


                         }


///////////////////////////////////////////////////////
                        ///////////////////////  Plot Exploded Hexes
                        ///////////////////////////////
/*
                        if (sc.stamp_centers)
                        {


                            if (u_dist > border_w && v_dist > border_w && w_dist > border_w)
                            {
                                acc_smoke_r.setValue(ijk, 0.95 );
                                acc_smoke_g.setValue(ijk, 0.95  );
                                acc_smoke_b.setValue(ijk, 0.95  );
                                acc_smoke_density.setValue(ijk, .2 );

                                acc_r.setValue(ijk, .75 );
                                acc_g.setValue(ijk, .75 );
                                acc_b.setValue(ijk, .75 );
                                acc_strength.setValue(ijk, .2);
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
                                acc_smoke_density.setValue(ijk, .01 );
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
                                acc_smoke_density.setValue(ijk, .01 );
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
                                acc_smoke_density.setValue(ijk, .01 );
                            }






                     ///   }


            //        }

                } */
            }
        }


        if (tetIsActive)
        {
            for (i = (int) curt_min[0]; i < (int) curt_max[0]; ++i) {
                for (j = (int) curt_min[1]; j < (int) curt_max[1]; ++j) {
                    for (k = (int) curt_min[2]; k < (int) curt_max[2]; ++k) {
                        Eigen::Vector3d p;
                        p << float(i), float(j), float(k);
                    
                    // Check if a voxel is in a tet ( they in same coordinate system now. )
                        Eigen::VectorXd param_val_tet;
                        bool pIsIn = pointInsideT(A, B, C, D, p, param_val_tet);
                        if ( pIsIn )
                        {
                            tetIsActive = true;
                            acc_r.setValue(ijk, param_val_tet(0) );
                            acc_g.setValue(ijk, param_val_tet(1) );
                            acc_b.setValue(ijk, param_val_tet(2) );
                                    // // acc_strength.setValue(ijk, unit_u);
                                    //     acc_strength.setValue(ijk, .5);

                                    //     acc_r.setValue(ijk, jac0 );
                                    //     acc_g.setValue(ijk, jac1 );
                                    //     acc_b.setValue(ijk, jac2 );
                                    // // acc_strength.setValue(ijk, unit_u);
                            acc_strength.setValue(ijk, .05);
                            acc_smoke_density.setValue(ijk, .1);
                            acc_smoke_r.setValue(ijk, .9 );
                            acc_smoke_g.setValue(ijk, .9 );
                            acc_smoke_b.setValue(ijk, .9 );
                        }
                    }
                }
            }
        }


    }


}

