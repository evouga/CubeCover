#ifndef VOXELUTILS_H
#define VOXELUTILS_H

#include <Eigen/Core>
#include <vector>
#include <glm/glm.hpp>
#include <openvdb/openvdb.h>
#include "ReadHexEx.h"



class sceneInfo
{
public:
    sceneInfo(const std::string hexexfile, const double sample_res );

    int nTets() const { return T.rows(); }
    // int nFaces() const { return F.rows(); }
    // int nEdges() const { return E.rows(); }


    Eigen::MatrixXd V;
    Eigen::MatrixXd V_pixel;
    Eigen::MatrixXi T;
    Eigen::MatrixXd param;
    Eigen::MatrixXd param_pixel;


    Eigen::Vector3d mesh_min;
    Eigen::Vector3d mesh_max;
    Eigen::Vector3d param_min;
    Eigen::Vector3d param_max;
};

bool pointInsideT(const Eigen::Vector3d& A, 
	              const Eigen::Vector3d& B, 
	              const Eigen::Vector3d& C, 
	              const Eigen::Vector3d& D, 
	              const Eigen::Vector3d& samplePoint, Eigen::VectorXd& textureCoordinate);


void stampParamView(openvdb::FloatGrid::Accessor acc_r, 
                    openvdb::FloatGrid::Accessor acc_g, 
                    openvdb::FloatGrid::Accessor acc_b, 
                    openvdb::FloatGrid::Accessor acc_smoke_r, 
                    openvdb::FloatGrid::Accessor acc_smoke_g, 
                    openvdb::FloatGrid::Accessor acc_smoke_b, 
                    openvdb::FloatGrid::Accessor acc_strength, // for now blackbody and emission are the same
                    openvdb::FloatGrid::Accessor acc_smoke_density,
                    openvdb::Vec3SGrid::Accessor acc_smoke_color);

    

/*

namespace CubeCover
{
    class TetMeshConnectivity;
};

void extractIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& param,
    Eigen::MatrixXd& P,
    Eigen::MatrixXi& E,
    Eigen::MatrixXd& P2,
    Eigen::MatrixXi& E2,
    std::vector<int> badverts);


void extractPoints(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& param,
    std::vector<glm::vec3>& points,
    std::vector<std::array<double, 3>>& colors );
*/

#endif
