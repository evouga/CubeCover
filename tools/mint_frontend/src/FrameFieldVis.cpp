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


            // std::cout << "sanity check should be 1: " << b1.norm() << std::endl;

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
    int startTetId,
    std::vector<Eigen::Vector2i>& tree_traversal,
    std::vector<Eigen::Vector4i>& tree_traversal_metadata
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
    tree_traversal.clear();
    tree_traversal_metadata.clear();
    // std::vector<Eigen::Vector2i> tree_traversal;

    // std::vector<std::vector<Eigen::Vector2i>> // source vert id, sink vert id, 

    int cur_tid = startTetId;

    // BASE CASE. 
    int init_edge_id = mesh.tetEdge(cur_tid,0);
    Eigen::Vector4i init_edge_metadata = Eigen::Vector4i(cur_tid, 0, init_edge_id, 1);

    Eigen::Vector2i init_edge = Eigen::Vector2i(mesh.edgeVertex(init_edge_id,0),mesh.edgeVertex(init_edge_id,1));

    tree_traversal.push_back(init_edge);
    tree_traversal_metadata.push_back(init_edge_metadata);
    vert_visited[mesh.edgeVertex(init_edge_id,0)] = 1;
    vert_visited[mesh.edgeVertex(init_edge_id,1)] = 1;

    tet_added_to_queue[cur_tid] += 1;

    tet_queue.push_back(cur_tid);
    // std::deque<int>::iterator it = tet_queue.begin();

    int iter = 0;
    while ( !tet_queue.empty() ) 
    // while( it != tet_queue.end() )
    {
        // cur_tid = *it;

        cur_tid = tet_queue.front();
        tet_queue.pop_front();

    // tryAddNeighborTetToQueue(cur_tid, tet_visited);
        for (int tetVert = 0; tetVert < 4; tetVert++)
        {
            int opp_tet_id = mesh.tetOppositeVertex(cur_tid, tetVert);

            int cur_face_id = mesh.tetFace(cur_tid, tetVert); // this might be wrong

            // std::cout << "opp_tet_id " << opp_tet_id << " face id " << cur_face_id << "blah " << field.faceAssignment(cur_face_id).isIdentity() << std::endl;


            if (opp_tet_id > -1 && field.faceAssignment(cur_face_id).isIdentity())
            {
                // If not in queue add it.
                std::cout << "got here " << std::endl;

                if (tet_added_to_queue[opp_tet_id] == 0)
                {
                    tet_added_to_queue[opp_tet_id] = 1;
                    tet_queue.push_back(opp_tet_id);
                }
            }

        }

        // it++;
        iter++;
        // std::cout << "boooooop" << iter << std::endl;


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
                int eid0 = mesh.edgeVertex(cur_edge_id,0);
                int eid1 = mesh.edgeVertex(cur_edge_id,1);
                
                if (v1_val > 0)
                {
                    int eid0 = mesh.edgeVertex(cur_edge_id,1);
                    int eid1 = mesh.edgeVertex(cur_edge_id,0);
                    orientation = -1;
                }
                    
                Eigen::Vector4i cur_edge_metadata = Eigen::Vector4i(cur_tid, tetEdge, init_edge_id, orientation);
                tree_traversal.push_back( Eigen::Vector2i(eid0,eid1) );
                tree_traversal_metadata.push_back( cur_edge_metadata );
                // tree_traversal.push_back( Eigen::Vector2i(cur_edge_id, orientation) );

                vert_visited[mesh.edgeVertex( cur_edge_id, 0 )] += 1;
                vert_visited[mesh.edgeVertex( cur_edge_id, 1 )] += 1;

            }

            
        }

    }

}




void integrateFieldOnEdges(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    const std::vector<Eigen::Vector2i>& tree_traversal,
    const std::vector<Eigen::Vector4i>& tree_traversal_metadata,
    double period,
    Eigen::MatrixXd& integratedVals
)
{


    int ntets = mesh.nTets();
    int nverts = V.rows(); 
    int nfaces = mesh.nFaces();
    int vpf = field.vectorsPerFrame();

    // integratedVals;
    integratedVals.resize(nverts, 2*vpf);
    // integratedVals

    for (int i = 0; i < nverts-1; i++)
    {
        // Eigen::MatrixXd cur_frame = frameVectors.at(i);
        Eigen::Vector3d cur_edge = V.row(tree_traversal.at(i)[1]) - 
                                       V.row(tree_traversal.at(i)[0]);

        int cur_tet_id = tree_traversal_metadata.at(i)[0];
        Eigen::VectorXd src_vals = integratedVals.row(tree_traversal.at(i)[0]);
        Eigen::VectorXd edge_diffs = src_vals*0;
        for (int j = 0; j < 2*vpf; j++)
        {
            edge_diffs[j] = frameVectors.at(j).row(cur_tet_id).dot(cur_edge);
        }

        Eigen::VectorXd sink_vals = src_vals + edge_diffs;
        if ( period > 0 )
        {
            sink_vals = sink_vals * 1./ period;
            for (int j = 0; j < 2*vpf; j++)
            {
                sink_vals[j] = sink_vals[j] - std::floor(sink_vals[j]);
                if (sink_vals[j] < 0.)
                {
                    sink_vals[j] = 1. - sink_vals[j];
                }
            }
            sink_vals = sink_vals * period;

        }

        std::cout << sink_vals << std::endl;

        integratedVals.row(tree_traversal.at(i)[1]) = sink_vals;



    }



}



void projectVertScalarsToTetFrames(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    const Eigen::MatrixXd& integratedVals,
    std::vector<Eigen::MatrixXd>& ret_frames
)
{ //     Eigen::MatrixXd& ret_frames,

    int ntets = mesh.nTets();
    int nverts = V.rows(); 
    int nfaces = mesh.nFaces();
    // int vpf = field.vectorsPerFrame();
    int inner_iters = frameVectors.size();

    ret_frames.resize(inner_iters);
    for (int i = 0; i < inner_iters; i++)
    {
        ret_frames[i].resize(ntets, 3);
    }


    for(int i = 0; i < ntets; i++)
    {
        for(int j = 0; j < inner_iters; j++)
        {
            int v0 = mesh.tetVertex(i, 0);
            int v1 = mesh.tetVertex(i, 1);
            int v2 = mesh.tetVertex(i, 2);
            int v3 = mesh.tetVertex(i, 3);

            Eigen::Vector3d e1 = V.row(v1) - V.row(v0);
            Eigen::Vector3d e2 = V.row(v2) - V.row(v0);
            Eigen::Vector3d e3 = V.row(v3) - V.row(v0);

            Eigen::Matrix3d m; 
            m.row(0) = e1; m.row(1) = e2; m.row(2) = e3;
            Eigen::Vector3d b; 
            b(0) = integratedVals(v1,j) - integratedVals(v0,j);
            b(1) = integratedVals(v2,j) - integratedVals(v0,j);
            b(2) = integratedVals(v3,j) - integratedVals(v0,j);

            Eigen::Vector3d f = m.ldlt().solve(b);

            ret_frames[j].row(i) = f;
        }

    }
}