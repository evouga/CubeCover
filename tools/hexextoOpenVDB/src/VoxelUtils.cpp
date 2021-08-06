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


    this->mesh_min = mesh_min;
    this->mesh_max = mesh_max;
    this->param_min = param_min;
    this->param_max = param_max;
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


void stampParamView(openvdb::FloatGrid::Accessor acc_r, 
                    openvdb::FloatGrid::Accessor acc_g, 
                    openvdb::FloatGrid::Accessor acc_b, 
                    openvdb::FloatGrid::Accessor acc_smoke_r, 
                    openvdb::FloatGrid::Accessor acc_smoke_g, 
                    openvdb::FloatGrid::Accessor acc_smoke_b, 
                    openvdb::FloatGrid::Accessor acc_strength, // for now blackbody and emission are the same
                    openvdb::FloatGrid::Accessor acc_smoke_density,
                    openvdb::Vec3SGrid::Accessor acc_smoke_color)
{




}