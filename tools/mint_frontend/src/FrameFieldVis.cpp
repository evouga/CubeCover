#include "FrameFieldVis.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <iostream>
#include <Eigen/Geometry> 
#include <deque>

void buildFrameVectors(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    double scale,
    Eigen::MatrixXd& centroids,
    std::vector<Eigen::MatrixXd>& frameVectors   
)
{
    int ntets = mesh.nTets();
    int vpf = field.vectorsPerFrame();
    centroids.resize(ntets, 3);
    frameVectors.resize(2*vpf);
    for (int i = 0; i < 2 * vpf; i++)
    {
        frameVectors[i].resize(ntets, 3);
    }

    for (int i = 0; i < ntets; i++)
    {
        Eigen::Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 4; j++)
        {
            int vert = mesh.tetVertex(i, j);
            centroid += V.row(vert).transpose();            
        }
        centroid /= 4.0;
        centroids.row(i) = centroid.transpose();

        for (int j = 0; j < vpf; j++)
        {
            Eigen::Vector3d vec = field.tetFrame(i).row(j).transpose();
            frameVectors[2 * j].row(i) = scale * vec.transpose();
            frameVectors[2 * j + 1].row(i) = -scale * vec.transpose();
        }
    }
}



void computePerVectorCurl(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    Eigen::MatrixXd& splitCurl
)
{
    int ntets = mesh.nTets();
    int nverts = V.rows(); 
    int nfaces = mesh.nFaces();
    int vpf = field.vectorsPerFrame();

    // for (int i = 0; i < ntets; i++)
    // { 


    // }

    splitCurl.resize(ntets, 6);

    // Eigen::MatrixXd& 

    for (int i = 0; i < nfaces; i++)
    { 
        if (!mesh.isBoundaryFace(i) && field.faceAssignment(i).isIdentity())
        {
            int tid0 = mesh.faceTet(i, 0);
            int tid1 = mesh.faceTet(i, 1);

            Eigen::Vector3d v0 = V.row(mesh.faceVertex(i, 0));
            Eigen::Vector3d v1 = V.row(mesh.faceVertex(i, 1));
            Eigen::Vector3d v2 = V.row(mesh.faceVertex(i, 2));
            
            Eigen::Vector3d fn = (v1 - v0).cross(v2 - v0);

            Eigen::Vector3d b0 = fn.cross(v1-v0);
            Eigen::Vector3d b1 = fn.cross(b0);
            b0.normalize();
            b1.normalize();


            std::cout << "sanity check should be 1" << b1.norm() << std::endl;

            for (int j = 0; j < 2*vpf; j++)
            {
                Eigen::Vector3d vf0 = frameVectors[j].row(tid0);
                Eigen::Vector3d vf1 = frameVectors[j].row(tid1);

                double faceCurl = (vf0.dot(b0) - vf1.dot(b0)) + (vf0.dot(b1) - vf1.dot(b1));
                splitCurl(tid0,j) += faceCurl*faceCurl;
                splitCurl(tid1,j) += faceCurl*faceCurl;

            }  
        }

        // for (int j = 0; j < vpf; j++)
        // {

        // }  

    }


   




}


/* 
This function takes in a source id, and does breadth first search on the tets and 
builds up an oriented list of edges this way.  


*/
void makeEdgeSpanningTree(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    int startTetId,
    Eigen::MatrixXd& splitCurl
)
{
    int ntets = mesh.nTets();
    int nverts = V.rows(); 
    int nfaces = mesh.nFaces();
    int vpf = field.vectorsPerFrame();

    std::vector<int> tet_added_to_queue(ntets, 0);
    std::vector<int> vert_visited(nverts, 0);

    std::deque<int> tet_queue;
    std::vector<Eigen::Vector2i> edge_queue; // intexed as tetid, tetedge 

    std::vector<Eigen::Vector2i> tree_traversal;

    // std::vector<std::vector<Eigen::Vector2i>> // source vert id, sink vert id, 

    int cur_tid = startTetId;

    // BASE CASE. 
    int init_edge_id = mesh.tetEdge(cur_tid,0);
    Eigen::Vector2i init_edge = Eigen::Vector2i(init_edge_id,1);
    tree_traversal.push_back(init_edge);
    vert_visited[mesh.edgeVertex(init_edge_id,0)] = 1;
    vert_visited[mesh.edgeVertex(init_edge_id,1)] = 1;

    tet_added_to_queue[cur_tid] += 1;

    tet_queue.push_back(cur_tid);


    while ( !tet_queue.empty() ) 
    {
        cur_tid = tet_queue.front();
        tet_queue.pop_front();

    // tryAddNeighborTetToQueue(cur_tid, tet_visited);
        for (int tetVert = 0; tetVert < 4; tetVert++)
        {
            int opp_tet_id = mesh.tetOppositeVertex(cur_tid, tetVert);

            int cur_face_id = mesh.tetFace(cur_tid, tetVert); // this might be wrong

            if (opp_tet_id > -1 && field.faceAssignment(cur_face_id).isIdentity())
            {
                // If not in queue add it.
                if (tet_added_to_queue[cur_tid] == 0)
                {
                    tet_added_to_queue[cur_tid] = 1;
                    tet_queue.push_back(cur_tid);
                }
            }

        }



    // Process current tet
    // i.e. add edges to tree 
    // and do other accounting 

    // tryAddEdge(tree_traversal, tet_visited, vert_visited,tet_queue);
        for (int tetEdge = 0; tetEdge < 6; tetEdge++)
        {
            int cur_edge_id = mesh.tetEdge(cur_tid, tetEdge);

            int v0_val = vert_visited[mesh.edgeVertex( cur_edge_id, 0 )];
            int v1_val = vert_visited[mesh.edgeVertex( cur_edge_id, 1 )];

            int orientation = 1;

            if ( ( v0_val > 0 && v1_val == 0) || (v0_val == 0 && v1_val > 0) )
            {
                if (v1_val > 0)
                    orientation = -1;

                tree_traversal.push_back( Eigen::Vector2i(cur_edge_id, orientation) );

                vert_visited[mesh.edgeVertex( cur_edge_id, 0 )] += 1;
                vert_visited[mesh.edgeVertex( cur_edge_id, 1 )] += 1;

            }

            
        }

    }

// finish increment current tet id 




    // for (int tetFace = 0; tetFace < 4; tetFace++)
    // {
    //     // tryAddNeighborTetToQueue(cur_tid, tet_visited);

    //     mesh.tetEdge

    //     for (int faceEdge = 0; faceEdge < 3; faceEdge++)
    //     {
    //         tryAddEdge(tree_traversal, tet_visited, vert_visited,tet_queue);
    //     }
    // }


    for (int i = 0; i<ntets; i++)
    {
        // add neighbors to queue 
        for (int j = 0; j < 4; j++)
        {


        }

    }


    // for (int i = 0; i < ntets; i++)
    // { 


    // }

    splitCurl.resize(ntets, 6);

    // Eigen::MatrixXd& 

    for (int i = 0; i < nfaces; i++)
    { 
        if (!mesh.isBoundaryFace(i) && field.faceAssignment(i).isIdentity())
        {
            int tid0 = mesh.faceTet(i, 0);
            int tid1 = mesh.faceTet(i, 1);

            Eigen::Vector3d v0 = V.row(mesh.faceVertex(i, 0));
            Eigen::Vector3d v1 = V.row(mesh.faceVertex(i, 1));
            Eigen::Vector3d v2 = V.row(mesh.faceVertex(i, 2));
            
            Eigen::Vector3d fn = (v1 - v0).cross(v2 - v0);

            Eigen::Vector3d b0 = fn.cross(v1-v0);
            Eigen::Vector3d b1 = fn.cross(b0);
            b0.normalize();
            b1.normalize();


            std::cout << "sanity check should be 1" << b1.norm() << std::endl;

            for (int j = 0; j < 2*vpf; j++)
            {
                Eigen::Vector3d vf0 = frameVectors[j].row(tid0);
                Eigen::Vector3d vf1 = frameVectors[j].row(tid1);

                double faceCurl = (vf0.dot(b0) - vf1.dot(b0)) + (vf0.dot(b1) - vf1.dot(b1));
                splitCurl(tid0,j) += faceCurl*faceCurl;
                splitCurl(tid1,j) += faceCurl*faceCurl;

            }  
        }

        // for (int j = 0; j < vpf; j++)
        // {

        // }  

    }


   




}









