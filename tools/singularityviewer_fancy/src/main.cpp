#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ReadFrameField.h"
#include <iostream>
#include <fstream>
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include "SingularCurveNetwork.h"
#include "readMeshFixed.h"
#include "polyscope/surface_mesh.h"
#include "FrameFieldVis.h"
#include "polyscope/point_cloud.h"
#include <random>
#include "WriteFrameField.h"
#include "ReadFancyData.h"

#include "polyscope/volume_mesh.h"
// #include "TetMeshConnectivity.h"
// #include "FrameField.h"
#include "CurlCorrection.h"



int main(int argc, char *argv[])
{
    // if (argc != 3)
    // {
    //     std::cerr << "Usage: singularityviewer_fancy (file_path_no_extension)" << std::endl;
    //     return -1;
    // }

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    
    // std::string slug = argv[1];

    // std::string meshfile = slug + ".mesh";
    // std::string frafile = slug + ".fra";

    std::string meshfile = "disk_3480_tets_gl3.mesh";
    std::string frafile = "disk_3480_tets_gl3.fra";


    bool showViz = true;
    std::string badverts;

    if (argc >= 4)
    {
        badverts = argv[3];
        if (badverts != "1")
            showViz = false;
    }

    std::cout << badverts << std::endl;


    Eigen::MatrixXi F;
    if (!CubeCover::readMESH(meshfile, V, T, F))
        return -1;

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    if (!CubeCover::readFrameField(frafile, "", T, frames, assignments, true))
        return -1;

    CubeCover::TetMeshConnectivity mesh(T);
    
    CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return -1;



    std::cout << "No face assignments provided, recomputing: ";
    std::cout.flush();
    field->computeLocalAssignments();
    std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;

    
    field->combAssignments();

    Eigen::MatrixXd Pblack;
    Eigen::MatrixXi Eblack;
    Eigen::MatrixXd Pblue;
    Eigen::MatrixXi Eblue;
    Eigen::MatrixXd Pgreen;
    Eigen::MatrixXi Egreen;
    extractSingularCurveNetwork(V, mesh, *field, Pgreen, Egreen, Pblue, Eblue, Pblack, Eblack);

    Eigen::MatrixXd centroids;
    std::vector<Eigen::MatrixXd> framefieldvecs;
    buildFrameVectors(V, mesh, *field, 1.0, centroids, framefieldvecs);

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

    // visualize the seams

    std::vector<int> seamfaces;

    int ninverted = 0;
    int nnontrivial = 0;

    for (int i = 0; i < nfaces; i++)
    {
        if (!field->faceAssignment(i).isIdentity())
        {
            seamfaces.push_back(i);
            nnontrivial++;
        }

        if (field->faceAssignment(i).orientation() == -1)
            ninverted++;
    }

    std::cout << "Non-identity face assignments: " << nnontrivial << std::endl;
    if (ninverted > 0)
    {
        std::cout << "Warning: " << ninverted << " face assignments are orientation-reversing" << std::endl;
    }

    int nseamtris = seamfaces.size();

    Eigen::MatrixXd seamV(3 * nseamtris, 3);
    Eigen::MatrixXi seamF(nseamtris, 3);
    for (int i = 0; i < nseamtris; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            seamF(i, j) = 3 * i + j;
            seamV.row(3 * i + j) = V.row(mesh.faceVertex(seamfaces[i], j));
        }
    }


///////////////////
///  New stuff
///////////////////
// 
//     std::string edges_mintfile = slug + "_edges.txt";
//     std::string edgecurlfile = slug + "_perEdgeCurl.txt";

//     Eigen::MatrixXi edgesdata_mint; 
//     Eigen::MatrixXd edgeCurl; 
//     CubeCover::readEdgeCurl(edges_mintfile, edgecurlfile, edgesdata_mint, edgeCurl);



    int vpf = field->vectorsPerFrame();


    Eigen::SparseMatrix<double> C;
    // buildCurlMatrix(field->vectorsPerFrame(), V, field->meshConnectivity(), C);
// URG throwing this here for now, weird linking issue
    // CubeCover::TetMeshConnectivity mesh = field->meshConnectivity();

//     {

//         std::vector<Eigen::Triplet<double> > Ccoeffs;

//         int nfaces = mesh.nFaces();
//         int row = 0;
//         for (int i = 0; i < nfaces; i++)
//         {
//             int t0 = mesh.faceTet(i, 0);
//             int t1 = mesh.faceTet(i, 1);
//             if (t0 == -1 || t1 == -1)
//                 continue;
//             int vertids[3];
//             for (int j = 0; j < 3; j++)
//             {
//                 vertids[j] = mesh.faceVertex(i, j);
//             }

//             Eigen::Vector3d pts[3];
//             for (int j = 0; j < 3; j++)
//             {
//                 pts[j] = V.row(vertids[j]).transpose();
//             }

//             Eigen::Vector3d n = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
//             n.normalize();

//             Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
//             for (int j = 0; j < vpf; j++)
//             {
//                 for (int k = 0; k < 3; k++)
//                 {
//                     for (int l = 0; l < 3; l++)
//                     {
//                         Ccoeffs.push_back({ row + 3 * j + k, 3 * vpf * t0 + 3 * j + l, P(k,l) });
//                         Ccoeffs.push_back({ row + 3 * j + k, 3 * vpf * t1 + 3 * j + l, -P(k,l) });
//                     }
//                 }
//             }
//             row += 3 * vpf;
//         }

    
// // std::cout << row << std::endl;
//         C.resize(row, 3 * vpf * mesh.nTets());
//         C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
    
//     }

    buildCurlMatrix(vpf,  V, mesh, C);
    

    int ntets = field->meshConnectivity().nTets();

    Eigen::VectorXd unrolled(ntets * vpf * 3);

    double minval = 100900000;
    double maxval = -100000000;
    for (int i = 0; i < ntets; i++)
    {
        for (int j = 0; j < vpf; j++)
        {
            unrolled.segment<3>(3 * vpf * i + 3 * j) = field->tetFrame(i).col(j);

            for (int k = 0; k < 3 ; k++)
            {   
                auto blah = field->tetFrame(i).col(j);
                if ( abs(blah(k)) > maxval)
                    maxval = abs(blah(k));
                if (abs(blah(k)) < minval)
                    minval = abs(blah(k));
            }

            

        }
    }

    std::cout << "min val" << minval << std::endl;
    std::cout << "max val" << maxval << std::endl;


    std::cout << unrolled.size() << std::endl;
    std::cout << "rows " <<C.rows() << "cols " << C.cols() << std::endl;
    std::cout << mesh.nTets() << std::endl;
    std::cout << vpf << std::endl;
    std::cout << "T.rows" << T.rows() <<std::endl;


    Eigen::VectorXd all_curl = C * unrolled;

    Eigen::VectorXd per_tet_sum_curl = Eigen::VectorXd::Zero(ntets);

    // std::cout << all_curl << std::endl;

    int row = 0;
    for (int i = 0; i < nfaces; i++)
    {
        int t0 = mesh.faceTet(i, 0);
        int t1 = mesh.faceTet(i, 1);
        if (t0 == -1 || t1 == -1)
            continue;

        for (int j = 0; j < vpf; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                    int blockId = row;
                    per_tet_sum_curl(t0) += abs(all_curl(blockId + 3*j + k));
                    per_tet_sum_curl(t1) += abs(all_curl(blockId + 3*j + k));
            }
        }
        row += 3 * vpf;
    }


std::cout << "allCurl.max" << all_curl.maxCoeff()  <<std::endl;
std::cout << "allCurl.max" <<  all_curl.minCoeff() <<std::endl;
std::cout << "unrolled.max" << unrolled.maxCoeff()  <<std::endl;
std::cout << "unrolled.max" <<  unrolled.minCoeff() <<std::endl;

    double unrolled_min = 1000000;
    for(int i = 0; i < unrolled.size(); i++ )
    {
        if ( unrolled_min > abs(unrolled(i)))
            unrolled_min = abs(unrolled(i));

    }

std::cout << "unrolled abs min" << unrolled_min << std::endl;

    double all_curl_min = 1000000;
    for(int i = 0; i < all_curl.size(); i++ )
    {
        if ( all_curl_min > abs(all_curl(i)))
            all_curl_min = abs(all_curl(i));

        if ( abs(all_curl(i)) < 1e-13 && abs(all_curl(i)) > 0)
            std::cout << "FISHY " << i << " " << abs(all_curl(i)) << std::endl;

    }

std::cout << "all_curl abs min" << all_curl_min << std::endl;


    double per_tet_sum_curl_min = 1000000;
    for(int i = 0; i < per_tet_sum_curl.size(); i++ )
    {
        if ( per_tet_sum_curl_min > abs(per_tet_sum_curl(i)))
            per_tet_sum_curl_min = abs(per_tet_sum_curl(i));

        if ( abs(per_tet_sum_curl(i)) < 1e-10 && abs(per_tet_sum_curl(i)) > 0)
            std::cout << "FISHY " << i << " " << abs(per_tet_sum_curl(i)) << std::endl;

    }

std::cout << "per_tet_sum abs min" << per_tet_sum_curl_min << std::endl;

    std::cout << "min frame entry: " << frames.minCoeff() << std::endl;
    std::cout << "max frame entry: " << frames.maxCoeff() << std::endl;


    
// int ntets = field->meshConnectivity().nTets();
// Eigen::VectorXd per_tet_sum_curl = Eigen::VectorXd::Zero(ntets);


///////////////////
///  Polyscope
///////////////////
// 


    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(0.0, 1.0);


        polyscope::init();

        auto *tetc = polyscope::registerPointCloud("Centroids", centroids);
        glm::vec3 dotcolor(0.1, 0.1, 0.1);
        tetc->setPointColor(dotcolor);
        tetc->setPointRadius(0.001);


        int vpf_polyscope = framefieldvecs.size();
        for (int i = 0; i < vpf_polyscope; i++)
        {
            std::stringstream ss;
            ss << "Frame Vector " << i;
            auto *vf = tetc->addVectorQuantity(ss.str(), framefieldvecs[i]);
            vf->setVectorColor({ dist(rng),dist(rng),dist(rng) });
            vf->setVectorRadius(0.001);
            vf->setEnabled(true);
        }

        // auto *green = polyscope::registerCurveNetwork("Singular Curves (+1/4)", Pgreen, Egreen);
        // green->setColor({ 0.0,1.0,0.0 });

        // auto *blue = polyscope::registerCurveNetwork("Singular Curves (-1/4)", Pblue, Eblue);
        // blue->setColor({ 0.0,0.0,1.0 });

        // auto *black = polyscope::registerCurveNetwork("Singular Curves (irregular)", Pblack, Eblack);
        // black->setColor({ 0.0,0.0,0.0 });

        polyscope::registerTetMesh("tet_mesh", V, T);
        auto scalarQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("real_curl_sum", per_tet_sum_curl);
        scalarQ->setEnabled(true);

        auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
        psMesh->setTransparency(0.2);
        psMesh->setSurfaceColor({ 0.5,0.5,0.0 });

        auto* seammesh = polyscope::registerSurfaceMesh("Seam", seamV, seamF);
        seammesh->setSurfaceColor({ 0.0, 0.0, 0.0 });

        // visualize!
        polyscope::show();
    

    delete field;
}


    // int vpf = framefieldvecs.size();


//     Eigen::SparseMatrix<double> C;
//     // buildCurlMatrix(field->vectorsPerFrame(), V, field->meshConnectivity(), C);
// // URG throwing this here for now, weird linking issue
//     // CubeCover::TetMeshConnectivity mesh = field->meshConnectivity();

//     {

//         std::vector<Eigen::Triplet<double> > Ccoeffs;

//         int nfaces = mesh.nFaces();
//         int row = 0;
//         for (int i = 0; i < nfaces; i++)
//         {
//             int t0 = mesh.faceTet(i, 0);
//             int t1 = mesh.faceTet(i, 1);
//             if (t0 == -1 || t1 == -1)
//                 continue;
//             int vertids[3];
//             for (int j = 0; j < 3; j++)
//             {
//                 vertids[j] = mesh.faceVertex(i, j);
//             }

//             Eigen::Vector3d pts[3];
//             for (int j = 0; j < 3; j++)
//             {
//                 pts[j] = V.row(vertids[j]).transpose();
//             }

//             Eigen::Vector3d n = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
//             n.normalize();

//             Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
//             for (int j = 0; j < vpf; j++)
//             {
//                 for (int k = 0; k < 3; k++)
//                 {
//                     for (int l = 0; l < 3; l++)
//                     {
//                         Ccoeffs.push_back({ row + 3 * j + k, 3 * vpf * t0 + 3 * j + l, P(k,l) });
//                         Ccoeffs.push_back({ row + 3 * j + k, 3 * vpf * t1 + 3 * j + l, -P(k,l) });
//                     }
//                 }
//             }
//             row += 3 * vpf;
//         }

    
// // std::cout << row << std::endl;
//         C.resize(row, 3 * vpf * mesh.nTets());
//         C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
    
//     }



//     int ntets = field->meshConnectivity().nTets();

//     Eigen::VectorXd unrolled(ntets * vpf * 3);
//     for (int i = 0; i < ntets; i++)
//     {
//         for (int j = 0; j < vpf; j++)
//         {
//             unrolled.segment<3>(3 * vpf * i + 3 * j) = field->tetFrame(i).col(j);
//         }
//     }

//     std::cout << unrolled.size() << std::endl;
//     std::cout << C.size() << std::endl;
//     std::cout << mesh.nTets() << std::endl;
//     std::cout << vpf << std::endl;


//     Eigen::VectorXd all_curl = C.transpose() * unrolled;

//     Eigen::VectorXd per_tet_sum_curl = Eigen::VectorXd::Zero(ntets);


//     int row = 0;
//     for (int i = 0; i < nfaces; i++)
//     {
//         int t0 = mesh.faceTet(i, 0);
//         int t1 = mesh.faceTet(i, 1);
//         if (t0 == -1 || t1 == -1)
//             continue;

//         for (int j = 0; j < vpf; j++)
//         {
//             for (int k = 0; k < 3; k++)
//             {
//                 for (int l = 0; l < 3; l++)
//                 {
//                     per_tet_sum_curl(t0) += all_curl(3 * vpf * t0 + 3 * j + l);
//                     per_tet_sum_curl(t1) += all_curl(3 * vpf * t1 + 3 * j + l);
//                 }
//             }
//         }
//         row += 3 * vpf;
//     }

