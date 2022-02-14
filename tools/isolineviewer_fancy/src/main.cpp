#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <iostream>
#include "TetMeshConnectivity.h"
#include "ReadHexEx.h"
#include "polyscope/surface_mesh.h"
#include "ExtractIsolines.h"
#include <Eigen/Dense>

#include "FrameField.h"
#include "TraceStreamlines.h"
#include "ReadFrameField.h"
#include "readMeshFixed.h"
// #include "FrameFieldVis.h"
#include "WriteFrameField.h"


#include "polyscope/point_cloud.h"
#include "polyscope/volume_mesh.h"

int main(int argc, char *argv[])
{
    // if (argc != 2)
    // {
    //     std::cerr << "Usage: isolineviewer (.hexex file)" << std::endl;
    //     return -1;
    // }

    if (argc != 4)
    {
        std::cerr << "Usage: isolineviewer (.mesh file) (.fra) (.hexex file)" << std::endl;
        return -1;
    }


    std::string meshfile = argv[1];
    std::string frafile = argv[2];
    std::string hexexfile = argv[3];
    std::string permfile; 
    std::string smoothness_root = meshfile.substr(0, meshfile.size()-5);

    // std::cout << smoothness_root << std::endl;


     ///////////////////////    ///////////////////////
     /////////     Load frame field       /////////////
     ///////////////////////     //////////////////////


    Eigen::MatrixXi F;
    Eigen::MatrixXd V_mesh;
    Eigen::MatrixXi T_mesh;
    if (!CubeCover::readMESH(meshfile, V_mesh, T_mesh, F))
        return -1;

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    if (!CubeCover::readFrameField(frafile, permfile, T_mesh, frames, assignments, true))
        return -1;

    CubeCover::TetMeshConnectivity mesh(T_mesh);
    
    CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return -1;


    field->computeLocalAssignments();
    field->combAssignments();


     ///////////////////////    ///////////////////////
     /////////     Build Streamlines       /////////////
     ///////////////////////     //////////////////////
    std::vector<Streamline> traces;
    int max_iter_per_trace = 700;

    Eigen::VectorXi init_tet_ids;
    // init_tet_ids.resize(6);
    // init_tet_ids << 1, 200, 300, 400, 500, 600, 700, 800, 900, 1000;

    traceStreamlines(V_mesh, mesh, *field, init_tet_ids, max_iter_per_trace, traces);


     ///////////////////////    ///////////////////////
     /////////     Isolines stuff       /////////////
     ///////////////////////     //////////////////////
    ///////////////////////

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    
    Eigen::MatrixXd values;
    if (!CubeCover::readHexEx(hexexfile, V, T, values))
    {
        std::cerr << "error reading the .hexex file" << std::endl;
        return -1;
    }

    CubeCover::TetMeshConnectivity mesh_hexex(T);
    
    Eigen::MatrixXd P;
    Eigen::MatrixXi E;

    Eigen::MatrixXd P2;
    Eigen::MatrixXi E2;

    extractIsolines(V, mesh_hexex, values, P, E, P2, E2);

     ///////////////////////    ///////////////////////
     /////////     Boundary Mesh       /////////////
     ///////////////////////     //////////////////////
    
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

    polyscope::init();


    auto *psCurves = polyscope::registerCurveNetwork("Isolines", P, E);
    psCurves->setRadius(0.002);
    auto *psCurves2 = polyscope::registerCurveNetwork("Bad Isolines", P2, E2);
    psCurves2->setRadius(0.003);
    // auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
    // psMesh->setTransparency(0.2);


    ///////////////////////    ///////////////////////
    /////////     Show Streamlines       /////////////
    ///////////////////////     //////////////////////
    // std::vector<Streamline> traces;
    // int max_iter_per_trace = 100;
    // Eigen::VectorXi init_tet_ids;
    // init_tet_ids.resize(10);
    // init_tet_ids << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    // traceStreamlines(V_mesh, mesh, *field, init_tet_ids, max_iter_per_trace, traces);

    Eigen::VectorXd tet_colors;
    tet_colors.resize(T_mesh.rows());
    tet_colors.setZero();
    int ntets_mesh = T_mesh.rows();
    int ntraces = traces.size();

    std::vector<Eigen::Vector4i> streamline_tets;


    for (int tid = 0; tid < ntraces; tid++)
    {
 //       streamline_tets.clear();
        int nsteps = traces.at(tid).tetIds.size();
        std::cout << "traceId: " << tid << "nsteps: " << nsteps << std::endl;
        for (int i = 0; i < nsteps; i++ )
        {
            auto cur_trace = traces.at(tid);
            int cur_tet_id = cur_trace.tetIds.at(i);
            streamline_tets.push_back(T_mesh.row(cur_tet_id));
            // std::cout << "cur_tet_id: " << cur_tet_id << std::endl;
            tet_colors(cur_tet_id) = 1.;
        }

        // auto *streamline_tets_mesh = polyscope::registerTetMesh("frame_mesh" + std::to_string(tid), V_mesh, streamline_tets);
        // streamline_tets_mesh->setTransparency(.3);
        
    }

    Eigen::MatrixXd cur_points; 
    Eigen::MatrixXd points; 
    int nfam = 6;

    for (int cur_fam_id = 0; cur_fam_id < nfam; cur_fam_id++)
    {
        int iter = 0;
        std::vector<Eigen::Vector2d> cur_line_edges;
        std::vector<Eigen::Vector3d> cur_points;

        for (int tid = 0; tid < ntraces; tid++)
        {
            int tid_fam_id = tid % nfam;
            if (tid_fam_id == cur_fam_id)
            {
     //       streamline_tets.clear();
                int nsteps = traces.at(tid).points.size();
                
                Eigen::Vector3d first_edge = traces.at(tid).points.at(1) - traces.at(tid).points.at(0);
                double fen = first_edge.norm();

                bool addLast = true;
                for (int i = 0; i < nsteps; i++ )
                {
                    // Eigen::Vector3d edge = traces.at(tid).points.at(i) - traces.at(tid).points.at(i+1); 
                    // if (edge.norm() > fen * 5 )
                    // {
                    //     addLast = false;
                    //     break;
                    // }
                    cur_points.push_back( traces.at(tid).points.at(i) );

                }
                // if (addLast)
                // {
                //     cur_points.push_back( traces.at(tid).points.at(nsteps-1) );
                // }

                

                
                for (int i = 0; i < cur_points.size()-1; i++ )
                {
                    cur_line_edges.push_back(Eigen::Vector2d(iter, iter+1) );
                    iter++;
                    // Eigen::Vector3d edge = traces.at(tid).points.at(i) - traces.at(tid).points.at(i+1); 
                    // if ( edge.norm() > fen * 2)
                    // {
                    //     i = i + nsteps;
                    // }
                    // else
                    // {
                    //     cur_line_edges.push_back(Eigen::Vector2d(iter, iter+1) );
                    //     iter++;
                    // }
                  
                }
                iter++;
                // std::cout << "traceId: " << tid << "nsteps: " << nsteps << std::endl;
                // for (int i = 0; i < nsteps-1; i++ )
                // {
                //     auto cur_trace = traces.at(tid);
                //     int cur_tet_id = cur_trace.tetIds.at(i);
                //     streamline_tets.push_back(T_mesh.row(cur_tet_id));
                //     // std::cout << "cur_tet_id: " << cur_tet_id << std::endl;
                //     tet_colors(cur_tet_id) = 1.;
                // }


            }


            
        }

        auto *single_streamline = polyscope::registerCurveNetwork("streamline" + std::to_string(cur_fam_id), cur_points, cur_line_edges);
        single_streamline->setTransparency(1);
        single_streamline->setRadius(0.003);


    }



 //    for (int tid = 0; tid < ntraces; tid++)
 //    {
 //        int fam_id = tid % 6;
 //        cur_iter = iter_streamline(fam_id);
 // //       streamline_tets.clear();
 //        int nsteps = traces.at(tid).points.size();
 //        std::vector<Eigen::Vector2d> cur_edge;
 //        for (int i = 0; i < nsteps-1; i++ )
 //        {
 //            cur_edge.push_back(Eigen::Vector2d(i, i+1) );
 //        }
 //        std::cout << "traceId: " << tid << "nsteps: " << nsteps << std::endl;
 //        for (int i = 0; i < nsteps-1; i++ )
 //        {
 //            auto cur_trace = traces.at(tid);
 //            int cur_tet_id = cur_trace.tetIds.at(i);
 //            streamline_tets.push_back(T_mesh.row(cur_tet_id));
 //            // std::cout << "cur_tet_id: " << cur_tet_id << std::endl;
 //            tet_colors(cur_tet_id) = 1.;
 //        }

 //        auto *single_streamline = polyscope::registerCurveNetwork("streamline" + std::to_string(tid), traces.at(tid).points, cur_edge);
 //        single_streamline->setTransparency(.7);
 //        single_streamline->setRadius(0.005);
        
 //    }

    // auto *single_streamline = polyscope::registerCurveNetwork("streamline" + std::to_string(tid), traces.at(tid).points, cur_edge);
    // single_streamline->setTransparency(.7);
    // psCurves2->setRadius(0.001);


    std::cout << "ntraces: " << ntraces << "tet_colors " << ntets_mesh << std::endl;



    std::string cur_smoothness = "_dirichlet_" + std::to_string(16);
    std::string dirichlet_path = smoothness_root + cur_smoothness + ".txt";

    std::ifstream ifs(dirichlet_path);
    Eigen::VectorXd channel_smoothness;
    channel_smoothness.resize(ntets_mesh);
    for (int i = 0; i < ntets_mesh; i++)
    {
        ifs >> channel_smoothness(i);
    }

    polyscope::options::transparencyRenderPasses = 32;

    auto blah_mesh = polyscope::registerTetMesh("tet_mesh", V_mesh, T_mesh);
    auto scalarQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity(cur_smoothness, channel_smoothness);
    scalarQ->setEnabled(true);

    blah_mesh->setTransparency(0.2);

    // for (int i = 0; i < ntets_mesh; i++)
    // {

    // }


    // auto *streamline_tets_mesh = polyscope::registerTetMesh("frame_mesh", V_mesh, streamline_tets);
    // // auto volScalarQ = streamline_tets_mesh->addCellScalarQuantity("active_tets", tet_colors);
    // // volScalarQ->setEnabled(true);
    // // volScalarQ->setMapRange({0,1.});
    // // volScalarQ->setColorMap("blues");
    // streamline_tets_mesh->setTransparency(.5);




    ///////////////////////


    

    // visualize!
    polyscope::show();
}



 //    for (int tid = 0; tid < ntraces; tid++)
 //    {
 // //       streamline_tets.clear();
 //        int nsteps = traces.at(tid).points.size();
 //        std::vector<Eigen::Vector2d> cur_edge;
 //        for (int i = 0; i < nsteps-1; i++ )
 //        {
 //            cur_edge.push_back(Eigen::Vector2d(i, i+1) );
 //        }
 //        std::cout << "traceId: " << tid << "nsteps: " << nsteps << std::endl;
 //        for (int i = 0; i < nsteps-1; i++ )
 //        {
 //            auto cur_trace = traces.at(tid);
 //            int cur_tet_id = cur_trace.tetIds.at(i);
 //            streamline_tets.push_back(T_mesh.row(cur_tet_id));
 //            // std::cout << "cur_tet_id: " << cur_tet_id << std::endl;
 //            tet_colors(cur_tet_id) = 1.;
 //        }

 //        auto *single_streamline = polyscope::registerCurveNetwork("streamline" + std::to_string(tid), traces.at(tid).points, cur_edge);
 //        single_streamline->setTransparency(.7);
 //        single_streamline->setRadius(0.005);
        
 //    }