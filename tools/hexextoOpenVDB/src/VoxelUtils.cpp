#include "VoxelUtils.h"
#include "TetMeshConnectivity.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"



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
