#include "VoxelUtils.h"
#include "TetMeshConnectivity.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

sceneInfo::sceneInfo( const std::string hexexfile, const double sample_res )
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXd param;
    // std::cout << hexexfile << std::endl;
    if (!CubeCover::readHexEx(hexexfile, V, T, param))
    {
        std::cerr << "error reading the .hexex file" << std::endl;
  //      return -1;
    }


    this->V = V;
    this->T = T;
    this->param = param;
   // Compute tet bounding boxes 
    int ntets = T.rows();
    int nverts = V.rows();

    // For tet mesh 
    double BIG_NUM = 100000000000000000.0;

    Eigen::Vector3d mesh_min(BIG_NUM,BIG_NUM,BIG_NUM);
    Eigen::Vector3d mesh_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);

    // For parameterization 

    Eigen::Vector3d param_min(BIG_NUM,BIG_NUM,BIG_NUM);
    Eigen::Vector3d param_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);

    for (int i = 0; i < ntets; i++)
    {
        int idx = 0;
        for (int j = 0; j < 4; j++)
        {
            Eigen::Vector3d cur_point = V.row( T(i,j) );
            for (int k = 0; k < 3; k++)
            {
                if( cur_point(k) < mesh_min(k) )
                    mesh_min(k) = cur_point(k);
                if( cur_point(k) > mesh_max(k) )
                    mesh_max(k) = cur_point(k);
            }

            Eigen::Vector3d param_point = param.row( T(i,j) );
            for (int k = 0; k < 3; k++)
            {
                if( param_point(k) < param_min(k) )
                    param_min(k) = param_point(k);
                if( param_point(k) > param_max(k) )
                    param_max(k) = param_point(k);
            }
            
        }

    }

    // Compute Scale factors


    // move min of bounding box to (0,0,0), rescale tet embedding to morph into sampling domain.  

    Eigen::Vector3d rangebound = mesh_max - mesh_min;
    double scale_fac = rangebound.maxCoeff();
    double pixeltoparam_scale = scale_fac / sample_res; 
    double topixel_scale = sample_res / scale_fac; 
    // set the max bounding box dimension in the param domain = 1, and multiply by sample res.
    // This places the parametric domain into the full sampling space

    // scratch.  
 //   Eigen::MatrixXi topixel_scale = (rangebound * ( scale_fac / sample_res)).asDiagonal();
  //  Eigen::MatrixXi pixeltoparam_scale = (rangebound * ( sample_res / scale_fac )).asDiagonal();

    Eigen::Vector3d param_rangebound = param_max - param_min;

    // Set the scale in the parameter domain so that largest dimension is 1.
    double param_scale_fac = double(param_rangebound.maxCoeff());


    // Move mesh into pixel space


        // V_pixel == V_pixel_space.  I.e. it embeds the parameterization in an sample_res X sample_res X sample_res cube.  
    Eigen::MatrixXd V_pixel = V;
    for (int i = 0; i < nverts; i++)
    {
        Eigen::Vector3d scale_verts = V_pixel.row(i);
        V_pixel.row(i) = (scale_verts - mesh_min) * topixel_scale;


        // rescale and shift parameterization so that the parameterization domain contains the desired number of grid cells along the max dimension.
        for (int j = 0; j < 4; j++)
        {
            Eigen::Vector3d tmp = ( param.row(4*i+j) );
            tmp = tmp - param_min;
            param.row(4*i+j) = tmp;
            param.row(4*i+j) *= ( ( 1. )  / param_scale_fac); 
        }

    }
    Eigen::MatrixXd param_pixel = param * sample_res;
    // V_pixel = param_pixel;
    this->V_pixel = V_pixel;
    this->param_pixel = param_pixel;
    this->param_unitcell = param;


    this->mesh_min = mesh_min;
    this->mesh_max = mesh_max;
    this->param_min = param_min;
    this->param_max = param_max;

    std::cout << "param_min: " << param_min << " param_max: " << param_max << std::endl; 

    std::cout << "param_min: " << param_min * ( ( 1. )  / param_scale_fac) << 
                 " param_max: " << param_max * ( ( 1. )  / param_scale_fac) << std::endl; 
}

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

    bool isIn = (tC.maxCoeff() < 1.001 && tC.minCoeff() > -0.001);

    Eigen::Vector3d normal = (B-A).cross(C-A);

    isIn = isIn && ( tC.sum() < 1.0001); // tet has l1 constraint.  

    textureCoordinate = Eigen::VectorXd(4);
    textureCoordinate = paramPoint; 
    // std::cout << textureCoordinate.transpose() << std::endl;
    return isIn;

}


void stampParamView(sceneInfo sc,
                    openvdb::FloatGrid::Accessor acc_r, 
                    openvdb::FloatGrid::Accessor acc_g, 
                    openvdb::FloatGrid::Accessor acc_b, 
                    openvdb::FloatGrid::Accessor acc_smoke_r, 
                    openvdb::FloatGrid::Accessor acc_smoke_g, 
                    openvdb::FloatGrid::Accessor acc_smoke_b, 
                    openvdb::FloatGrid::Accessor acc_strength, // for now blackbody and emission are the same
                    openvdb::FloatGrid::Accessor acc_smoke_density,
                    openvdb::Vec3SGrid::Accessor acc_smoke_color)
{
    double BIG_NUM = 100000000000000000.0;
    double line_w = .05;
    double cells = .75;

    sc.param_pixel = sc.param_pixel * 1.; // rescale factor / max # of cells in a row.
    Eigen::MatrixXd param_gridcell = sc.param_unitcell * cells;


    double border_w = .2;
    double facet_w = .01;
    double halo_w = .01;
    double cosmic_background = 0.0000;

    int ntets = sc.nTets();

    for (int t = 0; t < ntets; t++)
    {
        // compute tet bbox 
        Eigen::Vector3d curt_min(BIG_NUM,BIG_NUM,BIG_NUM);
        Eigen::Vector3d curt_max(-BIG_NUM,-BIG_NUM,-BIG_NUM);
        for (int v_idx = 0; v_idx < 4; v_idx++)
        {
            int rowId = sc.T(t, v_idx);
            // Eigen::Vector3d cur_point = sc.V_pixel.row( 4*t + v_idx );

            Eigen::Vector3d cur_point = sc.param_pixel.row( 4*t + v_idx );
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
        Eigen::Vector3d A = sc.param_pixel.row( 4*t + 0 );
        Eigen::Vector3d B = sc.param_pixel.row( 4*t + 1 );
        Eigen::Vector3d C = sc.param_pixel.row( 4*t + 2 );
        Eigen::Vector3d D = sc.param_pixel.row( 4*t + 3 );




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
                        // interp_param_to_pixel *= cells;


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



                        acc_smoke_r.setValue(ijk, 0.5 );
                        acc_smoke_g.setValue(ijk, 0.5  );
                        acc_smoke_b.setValue(ijk, 0.5  );
                        acc_smoke_density.setValue(ijk, 0.000 );

                        

                        // Draw the glowing grid.

/////////////////////////////////////////////////////////////////////////////////////
                        // TODO: make distance function of radius here - refactor.  
                        /////////////////////////////////////////////////////////

                        double u_dist = std::min(unit_u, fabs(1. - unit_u) );
                        double v_dist = std::min(unit_v, fabs(1. - unit_v) );
                        double w_dist = std::min(unit_w, fabs(1. - unit_w) );

                        if (u_dist > border_w && v_dist > border_w && w_dist > border_w)
                        {
                            acc_smoke_r.setValue(ijk, 0.7 );
                            acc_smoke_g.setValue(ijk, 0.7  );
                            acc_smoke_b.setValue(ijk, 0.7  );
                            acc_smoke_density.setValue(ijk, 0.02 );
                        }




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


                            if ( u_dist < halo_w)
                            {

                                acc_strength.setValue(ijk,  .5); // half emissive?
         
                                acc_smoke_density.setValue(ijk, .3 );
                                acc_r.setValue(ijk, 0. );
                                acc_g.setValue(ijk, .5 );
                                acc_b.setValue(ijk, .5 );


                                acc_smoke_r.setValue(ijk, 0. );
                                acc_smoke_g.setValue(ijk, .9 );
                                acc_smoke_b.setValue(ijk, .9 );
                                // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                            }
                            else if ( v_dist < halo_w )
                            {
                                acc_strength.setValue(ijk,  .5); // half emissive?
         
                                acc_smoke_density.setValue(ijk, .3 );
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, 0. );
                                acc_smoke_b.setValue(ijk, .9 );

                                acc_r.setValue(ijk, .5 );
                                acc_g.setValue(ijk, 0. );
                                acc_b.setValue(ijk, .5 );


                            }
                            else if ( w_dist < halo_w )
                            {
                                acc_strength.setValue(ijk,  .5); // half emissive?
         
                                acc_smoke_density.setValue(ijk, .3 );
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, .9 );
                                acc_smoke_b.setValue(ijk, 0. );



                                acc_r.setValue(ijk, .5 );
                                acc_g.setValue(ijk, .5 );
                                acc_b.setValue(ijk, 0. );

                            }



                            acc_strength.setValue(ijk, .2);
                            acc_smoke_density.setValue(ijk, .1);



                        }
                        else{
                            acc_smoke_density.setValue(ijk, .00);

                            acc_smoke_r.setValue(ijk, .9 );
                            acc_smoke_g.setValue(ijk, .9  );
                            acc_smoke_b.setValue(ijk, .9  );

                            // acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( 0.f ), 
                            //                                          float( 0.f ), 
                            //                                          float( 0.f )) );


                            if ( u_dist < facet_w)
                            {

                                acc_strength.setValue(ijk,  .5); // half emissive?
         
                                acc_smoke_density.setValue(ijk, .3 );
                                acc_r.setValue(ijk, 0. );
                                acc_g.setValue(ijk, .5 );
                                acc_b.setValue(ijk, .5 );


                                acc_smoke_r.setValue(ijk, 0. );
                                acc_smoke_g.setValue(ijk, .9 );
                                acc_smoke_b.setValue(ijk, .9 );
                                // acc_strength.setValue(ijk, float( interp_param_to_pixel(0) ));
                            }
                            else if ( v_dist < facet_w )
                            {
                                acc_strength.setValue(ijk,  .5); // half emissive?
         
                                acc_smoke_density.setValue(ijk, .3 );
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, 0. );
                                acc_smoke_b.setValue(ijk, .9 );

                                acc_r.setValue(ijk, .5 );
                                acc_g.setValue(ijk, 0. );
                                acc_b.setValue(ijk, .5 );


                            }
                            else if ( w_dist < facet_w )
                            {
                                acc_strength.setValue(ijk,  .5); // half emissive?
         
                                acc_smoke_density.setValue(ijk, .3 );
                                acc_smoke_r.setValue(ijk, .9 );
                                acc_smoke_g.setValue(ijk, .9 );
                                acc_smoke_b.setValue(ijk, 0. );



                                acc_r.setValue(ijk, .5 );
                                acc_g.setValue(ijk, .5 );
                                acc_b.setValue(ijk, 0. );

                            }



                            acc_strength.setValue(ijk, .2);
                            acc_smoke_density.setValue(ijk, .05);
                        }


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

                }
            }
        }
    }


}


void stampLatticeView(sceneInfo sc,
                    openvdb::FloatGrid::Accessor acc_r, 
                    openvdb::FloatGrid::Accessor acc_g, 
                    openvdb::FloatGrid::Accessor acc_b, 
                    openvdb::FloatGrid::Accessor acc_smoke_r, 
                    openvdb::FloatGrid::Accessor acc_smoke_g, 
                    openvdb::FloatGrid::Accessor acc_smoke_b, 
                    openvdb::FloatGrid::Accessor acc_strength, // for now blackbody and emission are the same
                    openvdb::FloatGrid::Accessor acc_smoke_density,
                    openvdb::Vec3SGrid::Accessor acc_smoke_color)
{

    double BIG_NUM = 100000000000000000.0;
    double line_w = .05;
    double border_w = .2;

    int ntets = sc.nTets();

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

                        // // This plots the param param per point in parameter space 
                        // acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( unit_v ) + float( unit_w ), 
                        //                                               float( unit_w ) + float( unit_u ), 
                        //                                               float( unit_v ) + float( unit_u )) );

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

                            // acc_smoke_color2.setValue(ijk, openvdb::Vec3s(float( 0.f ), 
                            //                                          float( 0.f ), 
                            //                                          float( 0.f )) );
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

}