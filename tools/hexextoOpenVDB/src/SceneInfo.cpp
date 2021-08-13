#include <vector>
#include <Eigen/Core>

#include "SceneInfo.h"
#include "ReadHexEx.h"
#include <iostream>
#include <fstream>

#include <openvdb/openvdb.h>

// namespace CubeCover {

SceneInfo::SceneInfo( const std::string hexexfile, const double sample_res )
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

    this->cells = 1;

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
 //   if ( V_curr == embedding::WORLD_SPACE )
    Eigen::Vector3d rangebound = mesh_max - mesh_min;



    double scale_fac = rangebound.maxCoeff();

    double topixel_scale = sample_res / scale_fac; 
    // set the max bounding box dimension in the param domain = 1, and multiply by sample res.
    // This places the parametric domain into the full sampling space



    Eigen::Vector3d param_rangebound = param_max - param_min;

    // Set the scale in the parameter domain so that largest dimension is 1.
    double param_scale_fac = double(param_rangebound.maxCoeff());


    // Move mesh into pixel space


        // V_pixel == V_pixel_space.  I.e. it embeds the parameterization in an sample_res X sample_res X sample_res cube.  
    Eigen::MatrixXd V_pixel = V;
    for (int i = 0; i < nverts; i++)
    {
        Eigen::Vector3d scale_verts = V_pixel.row(i);
        V_pixel.row(i) = (scale_verts - mesh_min) * (1./scale_fac);


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
    V_unitcell = V_pixel;
    V_pixel = V_unitcell * sample_res;
    // V_pixel = param_pixel;
    this->V_pixel = V_pixel;
    this->param_pixel = param_pixel;
    this->param_unitcell = param;


    this->mesh_min = mesh_min;
    this->mesh_max = mesh_max;
    this->param_min = param_min;
    this->param_max = param_max;

    this->stamp_grid = true;
    this->stamp_centers = true;

    std::cout << "param_min: " << param_min << " param_max: " << param_max << std::endl; 

    std::cout << "param_rescale: " << (param_max - param_min) * ( ( 1. )  / param_scale_fac) << std::endl; 
}



Eigen::Vector3d SceneInfo::get_V_pos(int t_idx, int v_idx)
{
    if( V_curr == embedding::PARAM_SPACE )
    {
        return param_pixel.row( 4*t_idx + v_idx );
    }
    else if( V_curr == embedding::WORLD_SPACE )
    {
        return V_pixel.row( T(t_idx, v_idx) );
    }
    return Eigen::Vector3d::Zero();

}



// }