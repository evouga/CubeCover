#ifndef SCENEINFO_H
#define SCENEINFO_H

#include <vector>
#include <Eigen/Core>

#include <openvdb/openvdb.h>

// namespace CubeCover {
enum embedding { PARAM_SPACE, WORLD_SPACE };

struct VDB_Grids {
    openvdb::FloatGrid::Ptr r;
    openvdb::FloatGrid::Ptr g; 
    openvdb::FloatGrid::Ptr b; 
    openvdb::FloatGrid::Ptr smoke_r; 
    openvdb::FloatGrid::Ptr smoke_g; 
    openvdb::FloatGrid::Ptr smoke_b; 
    openvdb::FloatGrid::Ptr strength; // for now blackbody and emission are the same
    openvdb::FloatGrid::Ptr smoke_density;
    openvdb::Vec3SGrid::Ptr smoke_color;

};

class SceneInfo
{
public:
    SceneInfo(const std::string hexexfile, const double sample_res );

    int nTets() const { return T.rows(); }
    int nVerts() const { return V.rows(); }


    int cells; // user provided parameter.  

    
    // int nFaces() const { return F.rows(); }
    // int nEdges() const { return E.rows(); }

    Eigen::Vector3d get_V_pos(int t_idx, int v_idx);

    embedding V_curr = PARAM_SPACE;
    bool stamp_grid;
    bool stamp_centers;
    VDB_Grids grids;

    Eigen::MatrixXi T;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXd V_unitcell;
    Eigen::MatrixXd V_pixel;
    // Eigen::MatrixXd V_curr; // should be in voxel coordinates

    Eigen::MatrixXd param;
    Eigen::MatrixXd param_pixel;
    Eigen::MatrixXd param_unitcell;


    Eigen::Vector3d mesh_min;
    Eigen::Vector3d mesh_max;
    Eigen::Vector3d param_min;
    Eigen::Vector3d param_max;


};
// };

#endif
