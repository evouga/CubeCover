#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ReadFrameField.h"
#include <iostream>
#include <fstream>
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include "readMeshFixed.h"
#include "FrameFieldVis.h"
#include <random>







int main(int argc, char *argv[])
{



	// todo add default out location
    if (argc != 5)
    {
        std::cerr << "Usage: proj_to_int_frames [int start idx] [.mesh] [.fra in filepath] [.fra out filepath]" << std::endl;
        return -1;
    }

    int start_node = std::stoi(argv[1]);
    std::string path_mesh = argv[2];
    std::string path_fra = argv[3];
    std::string path_outdir = argv[4];

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    Eigen::MatrixXi F;


    CubeCover::TetMeshConnectivity mesh;
    if (!CubeCover::readMESH(path_mesh, V, T, F))
    {
        return -1;
    }
    mesh = CubeCover::TetMeshConnectivity(T);

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    
    if (!CubeCover::readFrameField(path_fra, "", T, frames, assignments, true))
        return -1;
    
    // assignments.resize(0, 2 + 3);
    CubeCover::FrameField* field;

    field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return -1;

    // if (recomputeperms)
    // {
    std::cout << "No face assignments provided, recomputing: ";
    std::cout.flush();
    field->computeLocalAssignments();
    std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
    // }
    field->combAssignments();


    // extractSingularCurveNetwork(V, mesh, *field, Pgreen, Egreen, Pblue, Eblue, Pblack, Eblack);

    // tree_traversal.clear();
    Eigen::MatrixXd centroids;
    std::vector<Eigen::MatrixXd> framefieldvecs;
    Eigen::MatrixXd splitCurls;
    std::vector<Eigen::Vector2i> tree_traversal;
    std::vector<Eigen::Vector4i> tree_traversal_metadata;
    Eigen::MatrixXd treeIntegratedVals;
    std::vector<Eigen::MatrixXd> proj_framefieldvecs;


    buildFrameVectors(V, mesh, *field, 1.0, centroids, framefieldvecs);
    computePerVectorCurl(V, mesh, *field, framefieldvecs, splitCurls );
    makeEdgeSpanningTree(V, mesh, *field, start_node % mesh.nTets(), tree_traversal, tree_traversal_metadata);
    integrateFieldOnEdges(V, mesh, *field, framefieldvecs, tree_traversal, tree_traversal_metadata, 0, treeIntegratedVals);
    projectVertScalarsToTetFrames(V, mesh, *field, framefieldvecs,treeIntegratedVals, proj_framefieldvecs);




    {
        // write the frames
        std::ofstream ofs(path_outdir);
        if (!ofs)
        {
            return false;
        }
        ofs << "FRA 1" << std::endl;

        int nframes = field->meshConnectivity().nTets();
        int vpt = field->vectorsPerFrame();
        int type = 1;
        ofs << nframes << " " << vpt << " " << type << std::endl;

        for (int i = 0; i < nframes; i++)
        {
            for (int j = 0; j < vpt; j++)
            {
                Eigen::Vector3d v = proj_framefieldvecs[2*j].row(i);
                for (int k = 0; k < 3; k++)
                {
  
                    if (k != 0)
                        ofs << " ";
                    ofs << v(k);
                }
                ofs << std::endl;
            }
        }

        if (!ofs)
            return false;
    }


}


