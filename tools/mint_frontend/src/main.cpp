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


#if defined(WIN32) ||  defined(_WIN32)
	 gui->path_outdir = "C:\\Users\\fool\\Documents\\MATLAB\\integrable-moments\\output_frames_dir\\disk_idx_1_free_test_spectral_bound";

#elif defined(__unix__)
   gui->path_outdir = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/disk_idx_1_free_test_spectral_bound/";
   gui->path_outdir = "/home/josh/Documents/MATLAB/integrable-frames-3d/output_frames_dir/lap_circ_new_loc_proj_d_6_2/";
#endif

  std::cout << argc << std::endl;
    if (argc == 2)
    {
        std::cout << "Usage: pass in full path of folder you want to open, the gui is a bit bugged and I'm lazy to fix" << std::endl;
        std::string mint_out_dir = argv[1];
        std::cout << mint_out_dir.c_str() << std::endl;


        // thank you chatgpt for this code snippet 

        gui->path_outdir = new char[mint_out_dir.size() + 1]; 
        // Ensure null termination
        // gui->path_outdir[511] = '\0';

        // Check length and copy the string
        if (mint_out_dir.size() < 512) {
            std::strcpy(gui->path_outdir, mint_out_dir.c_str());
        } else {
            // Handle the case where the string is too long for the allocated buffer
            std::cout << "Error: mint_out_dir exceeds the allocated buffer size." << std::endl;
        }

        gui->path_outdir[mint_out_dir.size()] = '\0';




        std::cout << "bloop" << std::endl;
        // gui->path_outdir = mint_out_dir;
        // return -1;
    }
    if (argc > 2)
    {
        std::cerr << "Usage: call w no arguments to open default path or pass in mint output dir to visualize folder contents" << std::endl;
        return -1;
    }


 
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

