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

#include <experimental/filesystem> 

#include <igl/file_dialog_open.h>

#include "MintGUI.h"

// #include "../tinyfiledialogs/tinyfiledialogs.h"
//  https://stackoverflow.com/a/47651444   





std::string directory_path = "";


MintFrontend::MintGUI* gui;



int main(int argc, char *argv[])
{

    gui = new MintFrontend::MintGUI();

    
	/* 
	
	
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

    // std::string slug = argv[1];

    std::string meshfile = argv[1];
    std::string frafile = argv[2];
    std::string file_slug = argv[3];


    FileParts f = fileparts(file_slug);
    std::cout << file_slug<< "AAAAAAAAA " << f.path << "BBBBBBBBb " << f.name << " " << f.ext << " " << std::endl;

// Handle Mint_mesh as well, this is only for mint_fit
    std::string fit_target = f.path + "_analytic.fra";
    frafile = fit_target;

// fileparts(file_slug,fpath,fname,fext);


    // std::string meshfile = "disk_3480_tets_gl3.mesh";
    // std::string frafile = "disk_3480_tets_gl3.fra";


    // bool showViz = true;
    // std::string badverts;

    // if (argc >= 4)
    // {
    //     badverts = argv[3];
    //     if (badverts != "1")
    //         showViz = false;
    // }

    // std::cout << badverts << std::endl;


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

    // std::string perEdgeCurl_fid = file_slug + '_perEdgeCurl.txt';
    // std::string perEdgeCurl_fid = file_slug + '_interior_edges.txt';
    int ntets = field->meshConnectivity().nTets();


    std::string perFaceCurl_fid = file_slug + "_perFaceCurl.txt";
    std::string perFaceSmoothness_fid = file_slug + "_perFaceSmoothness.txt";
    std::string perFaceFrameMags_fid = file_slug + "_frameMags.txt";

    std::ifstream ifs(perFaceCurl_fid);
    Eigen::VectorXd curl_mint;
    curl_mint.resize(ntets);
    for (int i = 0; i < ntets; i++)
    {
        ifs >> curl_mint(i);
    }
    ifs.close();
    Eigen::VectorXd tmp = curl_mint;
    std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    std::pair<double,double> curl_mint_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));




    std::ifstream ifs2(perFaceSmoothness_fid);
    Eigen::VectorXd smoothness_mint;
    smoothness_mint.resize(ntets);
    for (int i = 0; i < ntets; i++)
    {
        ifs2 >> smoothness_mint(i);
    }
    ifs2.close();

    tmp = smoothness_mint;
    std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    std::pair<double,double> smoothness_mint_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));
    std::cout << "high smooth " << tmp(int( tmp.rows() * .9 ) ) << std::endl;


    std::ifstream ifs3(perFaceFrameMags_fid);
    Eigen::VectorXd framemags_mint;
    framemags_mint.resize(ntets);
    for (int i = 0; i < ntets; i++)
    {
        ifs3 >> framemags_mint(i);
    }
    ifs3.close();

    tmp = framemags_mint;
    std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    std::pair<double,double> framemags_mint_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));
    std::cout << "high smooth " << tmp(int( tmp.rows() * .9 ) ) << std::endl;



    // std::string small_fid = file_slug + "_smallframes.txt";
    // std::string mid_fid = file_slug + "_midframes.txt";
    // std::string big_fid = file_slug + "_largeframes.txt";

    // std::ifstream ifs4(small_fid);
    // Eigen::VectorXd small_frames;
    // small_frames.resize(ntets);
    // for (int i = 0; i < ntets; i++)
    // {
    //     ifs4 >> small_frames(i);
    // }
    // ifs4.close();

    // tmp = small_frames;
    // std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    // std::pair<double,double> smallframes_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));

    // std::ifstream ifs5(mid_fid);
    // Eigen::VectorXd mid_frames;
    // mid_frames.resize(ntets);
    // for (int i = 0; i < ntets; i++)
    // {
    //     ifs5 >> mid_frames(i);
    // }
    // ifs5.close();

    // tmp = mid_frames;
    // std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    // std::pair<double,double> midframes_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));

    // std::ifstream ifs6(big_fid);
    // Eigen::VectorXd big_frames;
    // small_frames.resize(ntets);
    // for (int i = 0; i < ntets; i++)
    // {
    //     ifs6 >> big_frames(i);
    // }
    // ifs6.close();

    // tmp = big_frames;
    // std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    // std::pair<double,double> bigframes_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));




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

    buildCurlMatrix(vpf,  V, *field, C);
    



    Eigen::VectorXd unrolled(ntets * vpf * 3);

    double minval = 100900000;
    double maxval = -100000000;
    for (int i = 0; i < ntets; i++)
    {
        for (int j = 0; j < vpf; j++)
        {
            unrolled.segment<3>(3 * vpf * i + 3 * j) = field->tetFrame(i).row(j);

            for (int k = 0; k < 3 ; k++)
            {   
                auto blah = field->tetFrame(i).row(j);
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

    tmp = per_tet_sum_curl;
    std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
    std::pair<double,double> curl_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));

    std::cout << "min frame entry: " << frames.minCoeff() << std::endl;
    std::cout << "max frame entry: " << frames.maxCoeff() << std::endl;


    
// int ntets = field->meshConnectivity().nTets();
// Eigen::VectorXd per_tet_sum_curl = Eigen::VectorXd::Zero(ntets);


///////////////////
///  Polyscope
///////////////////
// 

	*/

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    polyscope::options::automaticallyComputeSceneExtents = true;

    polyscope::options::autoscaleStructures = true;
    polyscope::options::autocenterStructures = true;
    // polyscope::state::lengthScale = 1.;
    // // polyscope::options::automaticallyComputeSceneExtents = true;
    // polyscope::state::boundingBox = 
    // std::tuple<glm::vec3, glm::vec3>{ {-1., -1., -1.}, {1., 1., 1.} };


    polyscope::init();
    polyscope::options::transparencyRenderPasses = 32;


    // gui->path_mesh = "/home/josh/Documents/MATLAB/integrable-frames-3d/meshes/TetWild/square_slab_1000.mesh";
    // gui->set_base_mesh();
    // gui->show_base_mesh();
  //  gui->path_outdir = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/disk_idx_1_free_test_spectral_bound/";

	gui->path_outdir = "C:\\Users\\fool\\Documents\\MATLAB\\integrable-moments\\output_frames_dir\\disk_idx_1_free_test_spectral_bound";
    gui->load_state_from_output_dir();

		/*

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

        auto *green = polyscope::registerCurveNetwork("Singular Curves (+1/4)", Pgreen, Egreen);
        green->setColor({ 0.0,1.0,0.0 });
        green->setTransparency(0.8);

        auto *blue = polyscope::registerCurveNetwork("Singular Curves (-1/4)", Pblue, Eblue);
        blue->setColor({ 0.0,0.0,1.0 });
        blue->setTransparency(0.8);

        auto *black = polyscope::registerCurveNetwork("Singular Curves (irregular)", Pblack, Eblack);
        black->setColor({ 0.0,0.0,0.0 });
        black->setTransparency(0.8);

        polyscope::registerTetMesh("tet_mesh", V, T);
        auto scalarQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("real_curl_sum", per_tet_sum_curl);
        scalarQ->setEnabled(true);
        scalarQ->setMapRange(curl_range);

        auto mintCurlQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("mint_curl_sum", curl_mint);
        mintCurlQ->setEnabled(false);
        mintCurlQ->setMapRange(curl_mint_range);

        auto mintSmoothQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("mint_smoothness_sum", smoothness_mint);
        mintSmoothQ->setEnabled(false);
        mintSmoothQ->setMapRange(smoothness_mint_range);

        auto mintFrameMagsQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("mint_framemags_sum", framemags_mint);
        mintFrameMagsQ->setEnabled(false);
        mintFrameMagsQ->setMapRange(framemags_mint_range);

        CubeCover::plot_tet_scalar_field(file_slug + "_perFaceCurl.txt", "tet_mesh", "FACE CURL");
        CubeCover::plot_tet_scalar_field(file_slug + "_smallframes.txt", "tet_mesh", "SMALL NORMS");
        CubeCover::plot_tet_scalar_field(file_slug + "_medframes.txt", "tet_mesh", "MEDIUM NORMS");
        CubeCover::plot_tet_scalar_field(file_slug + "_largeframes.txt", "tet_mesh", "BIGGEST NORMS");



        // auto smallFrameMagsQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("small_framemags", small_frames);
        // smallFrameMagsQ->setEnabled(false);
        // smallFrameMagsQ->setMapRange(smallframes_range);

        // auto midFrameMagsQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("mid_framemags", mid_frames);
        // midFrameMagsQ->setEnabled(false);
        // midFrameMagsQ->setMapRange(midframes_range);

        // auto bigFrameMagsQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity("big_framemags", big_frames);
        // bigFrameMagsQ->setEnabled(false);
        // bigFrameMagsQ->setMapRange(bigframes_range);

        

        polyscope::getVolumeMesh("tet_mesh")->setTransparency(0.2);


        auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
        psMesh->setTransparency(0.2);
        psMesh->setSurfaceColor({ 0.5,0.5,0.0 });
        psMesh->setEnabled(false);

        auto* seammesh = polyscope::registerSurfaceMesh("Seam", seamV, seamF);
        seammesh->setSurfaceColor({ 0.0, 0.0, 0.0 });
        seammesh->setEnabled(false);

		*/


		// polyscope::state::userCallback = myCallback;
        polyscope::state::userCallback = [](){gui->gui_callback();} ;
        
        // visualize!
        polyscope::show();
    

   // delete field;
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

