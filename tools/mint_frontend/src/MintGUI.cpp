#include <igl/file_dialog_open.h>
#include <igl/readOBJ.h>
#include <igl/bfs_orient.h>

#include "MintGUI.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

// #include <imgui/misc/cpp/imgui_stdlib.h>
// #include <polyscope/deps/imgui/imgui/imgui.h>
// #include <imgui/imgui.h>  
#include "TetMeshConnectivity.h"
#include "readMeshFixed.h"
#include "ReadMoments.h"

#include <misc/cpp/imgui_stdlib.h>
#include <algorithm>
#include <cctype>
#include <string>


#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;


// This is c++14 experimental feature.  In c++17 >= this is part of STD.
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;


namespace MintFrontend
{


static void HelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


	MintGUI::MintGUI()
	{
        path_mesh = new char[512];
        path_constraints = new char[512];
        path_outdir = new char[512];

        mesh = CubeCover::TetMeshConnectivity();
        moment_view_mode = Moments_To_Show::fourth;

        showBoundary = true;
        showInteriorTets = true;


        const char *args[] = {"x^4", "x^3 y", "x^3 z", "x^2 y^2", "x^2 y z", 
                              "x^2 z^2", "x y^3", "x y^2 z", "x y z^2", "x z^3", 
                              "y^4", "y^3 z", "y^2 z^2", "y z^3", "z^4"};
        std::vector<std::string> tmp(args, std::end(args));
        moment_labels = tmp;



        // exploded_spacing = 120.; // square

        exploded_spacing = 22.;

        sel_idx = -1;


        // Load config file here. 

        
	}

char * fileSelectSubroutine()
{
    std::string picked_file = igl::file_dialog_open();
    return strdup(picked_file.c_str());
    // return picked_file.data();
}






//! Using only text manipulation, splits a full path into component file parts
FileParts MintGUI::fileparts(const std::string &fullpath)
{
    using namespace std;

    size_t idxSlash = fullpath.rfind("/");
    if (idxSlash == string::npos) {
        idxSlash = fullpath.rfind("\\");
    }
    size_t idxDot = fullpath.rfind(".");

    FileParts fp;
    if (idxSlash != string::npos && idxDot != string::npos) {
        fp.path = fullpath.substr(0, idxSlash + 1);
        fp.name = fullpath.substr(idxSlash + 1, idxDot - idxSlash - 1);
        fp.ext  = fullpath.substr(idxDot);
    } else if (idxSlash == string::npos && idxDot == string::npos) {
        fp.name = fullpath;
    } else if (/* only */ idxSlash == string::npos) {
        fp.name = fullpath.substr(0, idxDot);
        fp.ext  = fullpath.substr(idxDot);
    } else { // only idxDot == string::npos
        fp.path = fullpath.substr(0, idxSlash + 1);
        fp.name = fullpath.substr(idxSlash + 1);
    }
    return fp;
}


void MintGUI::show_base_mesh()
{
    std::cout << "show_base_mesh" << std::endl;
    clear_polyscope_state();
    polyscope::options::automaticallyComputeSceneExtents = true;

    auto tet_mesh = polyscope::registerTetMesh("tet_mesh", V, T)->setEdgeWidth(0.5)->setTransparency(.7);
    auto surf_mesh = polyscope::registerSurfaceMesh("surf_mesh", V, bdryF)->setEdgeWidth(1)->setTransparency(.7);

    tet_mesh->rescaleToUnit();
    surf_mesh->rescaleToUnit();

    tet_mesh->resetTransform();
    surf_mesh->resetTransform();
    
        // std::cout << "V: " << V.size() << std::endl;
        // std::cout << "T: " << T << std::endl;

    polyscope::view::resetCameraToHomeView();
    //                 polyscope::state::boundingBox = 
    // std::tuple<glm::vec3, glm::vec3>{ {-1., -1., -1.}, {1., 1., 1.} };
}



void MintGUI::show_constraint_vals()
{

    std::cout << "show_constraint_vals" << std::endl;
    clear_polyscope_state();
    polyscope::options::automaticallyComputeSceneExtents = true;
    // polyscope::state::lengthScale = polyscope::state::lengthScale * 1/5.;

    for( int i = 0; i < 3; i++)
    {
        for( int j = 0; j < 5; j++)
        {
            int cur_id = i*5 + j;
            glm::vec3 shift(-j * exploded_spacing, -i * exploded_spacing, 0);

            if (showInteriorTets)
            {
                // auto cur_mesh = polyscope::registerTetMesh(moment_labels.at(cur_id) + "_tm_" + std::to_string(1000000 + cur_id), V, T);
                auto cur_mesh = polyscope::registerTetMesh(std::to_string(10000 + cur_id) + "____" +moment_labels.at(cur_id), V, T);


                cur_mesh->setEdgeWidth(0.5)->setTransparency(.7)->rescaleToUnit();
                cur_mesh->resetTransform();
                cur_mesh->translate(shift);
                Eigen::VectorXd tmp_mvals = M_curr.block(0,cur_id,mesh.nTets(),cur_id+1);
                cur_mesh->addCellScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
            }

            if (showBoundary)
            {
                // auto surf_mesh = polyscope::registerSurfaceMesh(moment_labels.at(cur_id) + "_sm_" + std::to_string(1000000 + cur_id), V, bdryF);
                auto surf_mesh = polyscope::registerSurfaceMesh(std::to_string(20000 + cur_id) + "____" + moment_labels.at(cur_id), V, bdryF);



                surf_mesh->setEdgeWidth(0.5)->setTransparency(.7)->rescaleToUnit();
                surf_mesh->resetTransform();
                surf_mesh->translate(shift);

                std::cout << M_curr.rows()-mesh.nTets() << " diff " << M_curr.rows()-mesh.nTets() - bdryF.rows() << std::endl;
                Eigen::VectorXd tmp_mvals = M_curr.block(mesh.nTets(),cur_id,bdryF.rows(),cur_id+1);
                surf_mesh->addFaceScalarQuantity("cur-moment", tmp_mvals)->setEnabled(true);
            }


            // std::cout << "tet_mesh_" + std::to_string(1000 + cur_id) + " " <<  glm::to_string(cur_mesh->getTransform()) << std::endl;
        }
    }
    // polyscope::registerTetMesh("tet_mesh", V, T)->setEdgeWidth(0.5)->setTransparency(.7);
    // polyscope::registerSurfaceMesh("surf_mesh", V, bdryF)->setEdgeWidth(1)->setTransparency(.7);
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = polyscope::state::lengthScale * 3.;

    polyscope::view::resetCameraToHomeView();
    // polyscope::state::boundingBox = 
    //     std::tuple<glm::vec3, glm::vec3>{ {-2.5, -1.5, -1.}, {2.5, 1.5, 1.} };
}








void MintGUI::clear_polyscope_state()
{

    std::cout << "clear_polyscope_state" << std::endl;
    polyscope::removeStructure("tet_mesh", false);
    polyscope::removeStructure("surf_mesh", false);


    for( int i = 0; i < 3; i++)
    {
        for( int j = 0; j < 5; j++)
        {
            int cur_id = i*5 + j;
            polyscope::removeStructure(std::to_string(10000 + cur_id) + "____" + moment_labels.at(cur_id), false);
        }
    }

        for( int i = 0; i < 3; i++)
    {
        for( int j = 0; j < 5; j++)
        {
            int cur_id = i*5 + j;
            polyscope::removeStructure(std::to_string(20000 + cur_id) + "____" + moment_labels.at(cur_id), false);
        }
    }

}









void MintGUI::set_base_mesh()
{
    if (!CubeCover::readMESH(path_mesh, V, T, F))
    {
        polyscope::warning("Unable to load selected mesh");
        mesh = CubeCover::TetMeshConnectivity();
    }

    else{
        mesh = CubeCover::TetMeshConnectivity(T);
        // make boundary mesh out of volume mesh
            // make a mesh out of all of the boundary faces

        auto parts = fileparts(path_mesh);
        
        Eigen::MatrixXd bdryV;
        Eigen::MatrixXi tmp;

        igl::readOBJ(parts.path + parts.name + ".obj", bdryV, bdryF);

        igl::bfs_orient(bdryF, bdryF, tmp);
        
        // int nbdry = 0;
        // int nfaces = mesh.nFaces();
        // for (int i = 0; i < nfaces; i++)
        // {
        //     if (mesh.isBoundaryFace(i))
        //         nbdry++;
        // }

        // bdryF.resize(nbdry, 3);
        // // std::cout << "nbdry" << nbdry << std::endl<< std::endl<< std::endl;
        // // std::cout << "bdryF " << bdryF.size() << std::endl;
        // int curidx = 0;
        // for (int i = 0; i < nfaces; i++)
        // {
        //     if (mesh.isBoundaryFace(i))
        //     {
        //         for (int j = 0; j < 3; j++)
        //         {
        //             bdryF(curidx, j) = mesh.faceVertex(i, j);
                    
        //         }
        //         // fix triangle orientations
        //         int tet = mesh.faceTet(i, 0);
        //         if (tet == -1)
        //         {
        //             std::swap(bdryF(curidx, 0), bdryF(curidx, 1));
        //         }
        //         curidx++;
        //     }
        // }

        
        show_base_mesh();
    }


}

void MintGUI::load_state_from_output_dir()
{
    folder_contents.clear();
    file_names.clear();
    for (const auto & entry : fs::directory_iterator(path_outdir))
    {
        // std::cout << entry.path() << std::endl;
        std::string tmp = entry.path();
        FileParts fp = fileparts(tmp);
        if (fp.ext == ".mom"){
            folder_contents.push_back( entry.path() );
            file_names.push_back( fp.name );
        }

        if (fp.ext == ".mesh")
        {
            const char* tmp_path_mesh = tmp.c_str();
            std::cout << tmp_path_mesh << std::endl;
            path_mesh = new char[512];
            strncpy(path_mesh, tmp_path_mesh, 512);
        }
            // folder_contents.push_back( entry.path() );
    }

    sel_idx = -1;

    std::sort(folder_contents.begin(), folder_contents.end());
    std::sort(file_names.begin(), file_names.end());

    // for (int i = 0; i < folder_contents.size(); i++)
    // {
    //     std::cout << folder_contents.at(i) << std::endl;

    // }

    set_base_mesh();
    show_base_mesh();

}




void MintGUI::gui_callback()
{
    	ImGui::PushItemWidth(1000); // Make ui elements 100 pixels wide,
							   // instead of full width. Must have 
							   // matching PopItemWidth() below.

//	ImGui::InputInt("num points", &nPts);             // set a int variable
//	ImGui::InputFloat("param value", &anotherParam);  // set a float variable

//	if (ImGui::Button("run subroutine")) {
		// executes when button is pressed
//		mySubroutine();
//	}
    // std::string blah;
    // std::string blah2;
    // ImGui::InputText("Mesh Path", blah);
    // ImGui::InputText("Boundary Constraints Path", blah2);



    // std::cout << rightWindowsWidth << std::endl;


    ///////////////////////////////////////////////////////////////////////////

    ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("1", "abs mesh path, or use picker button.", path_mesh, 512);   
    ImGui::PopItemWidth();
    ImGui::SameLine();
    HelpMarker("Mesh Path");
    ImGui::SameLine();


    


    if (ImGui::Button("Pick .mesh")) {
    // executes when button is pressed
        char* tmp_path_mesh = fileSelectSubroutine();
        path_mesh = new char[512];
        strncpy(path_mesh, tmp_path_mesh, 512);
        auto cur_path_parts = fileparts(path_mesh);
        std::cout << cur_path_parts.ext << std::endl;

        std::string data = cur_path_parts.ext;
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });


        std::cout << data << std::endl;
        if (data == ".mesh" )
        {
            set_base_mesh();
        }
        else 
        {
            polyscope::warning("Please pick a .mesh file to load.  Support for .obj coming eventually...");
        }
    }

    if (ImGui::Button("(re)load view mesh")) {
        auto cur_path_parts = fileparts(path_mesh);
        std::cout << cur_path_parts.ext << std::endl;

        std::string data = cur_path_parts.ext;
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });
        if (data == ".mesh" )
        {
            set_base_mesh();
        }
        else 
        {
            polyscope::warning("Please pick a .mesh file to load.  Support for .obj coming eventually...");
        }
    }



    ///////////////////////////////////////////////////////////////////////////

        ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("2", "enter rel or abs path, or use picker.", path_outdir, 512);
    ImGui::PopItemWidth();


    ImGui::SameLine();
    HelpMarker("Specify directory to load a previous run, or to choose output for next run");
    ImGui::SameLine();


    if (ImGui::Button("Pick mint output dir")) {
        char* tmp_path_outdir = fileSelectSubroutine();
        strncpy(path_outdir, tmp_path_outdir, 512);
    }

 
    
    if (ImGui::Button("(re)load or create directory")) {

      polyscope::warning("This directory did not exist, creating it");
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();


        //     char* tmp_path_mesh = fileSelectSubroutine();
        // strncpy(path_mesh, tmp_path_mesh, 512);
        // auto cur_path_parts = fileparts(path_mesh);
        // std::cout << cur_path_parts.ext << std::endl;

        // std::string data = cur_path_parts.ext;
        // std::transform(data.begin(), data.end(), data.begin(),
        //     [](unsigned char c){ return std::tolower(c); });


        // std::cout << data << std::endl;
        // if (data == ".mesh" )
        // {
        //     set_base_mesh();
        // }
        // else 
        // {
        //     polyscope::warning("Please pick a .mesh file to load.  Support for .obj coming eventually...");
        // }




    }

  
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);

    if (ImGui::TreeNode("Contents of chosen dir"))
    {
             
        ImGui::SetNextItemOpen(true, ImGuiCond_Once);

        if (ImGui::TreeNode("TODO: Select Mesh from directory"))
        {

            ImGui::TreePop();
        }
        
        ImGui::SetNextItemOpen(true, ImGuiCond_Once);

        if (ImGui::TreeNode("Select Moments to Visualize"))
        {
             

            for (int i = 0; i < file_names.size(); i++)
            {
                if (ImGui::Selectable(file_names.at(i).c_str(), sel_idx == i))
                {
                    sel_idx = i;
                    path_constraints = new char[512];
                    strncpy(path_constraints, folder_contents.at(i).c_str(), 512);
                    CubeCover::readMoments(path_constraints, M_curr, true);
                    // std::cout << path_constraints << std::endl;
                    show_constraint_vals();

                }

            }


            ImGui::PushItemWidth(300);
            ImGui::InputTextWithHint("3", "path to exploded moments", path_constraints, 512);
            ImGui::PopItemWidth();
            ImGui::SameLine();
            HelpMarker("Choose moments to visualize in a different directory");
            ImGui::SameLine();

            if (ImGui::Button("Pick .mom")) {
                char* tmp_path_constraints = fileSelectSubroutine();
                path_constraints = new char[512];
                strncpy(path_constraints, tmp_path_constraints, 512);
                show_constraint_vals();

                sel_idx = -1;

            }

            ImGui::TreePop();
        }
        ImGui::TreePop();

    }



    ///////////////////////////////////////////////////////////////////////////



    
    if (ImGui::Button("(re)Load Boundary")) {

      polyscope::warning("The chosen .bound file does not match the loaded mesh.");
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
      show_constraint_vals();
    }



    if (ImGui::Button("Compute and Load normal boundary")) {

        show_constraint_vals();

    }

    ImGui::SetNextItemOpen(true, ImGuiCond_Once);

    if (ImGui::TreeNode("Exploded GUI options"))
    {

        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Exploded Spacing", &exploded_spacing, 0.0f, 150.0f, "ratio = %.3f");
        ImGui::PopItemWidth();

        ImGui::Text("Which moments to show");
        if (ImGui::RadioButton("4th", moment_view_mode == Moments_To_Show::fourth))  
        { 
            moment_view_mode = Moments_To_Show::fourth; 
            // polyscope::state::lengthScale = 5.;
        } ImGui::SameLine();
        if (ImGui::RadioButton("2nd", moment_view_mode == Moments_To_Show::second))  { moment_view_mode = Moments_To_Show::second; } ImGui::SameLine();
        if (ImGui::RadioButton("both", moment_view_mode == Moments_To_Show::both))  { moment_view_mode = Moments_To_Show::both; } 



        if (ImGui::Checkbox("show boundary", &showBoundary)) { show_constraint_vals(); }ImGui::SameLine();
        if (ImGui::Checkbox("show interior tets", &showInteriorTets)) { show_constraint_vals(); }


        ImGui::TreePop();

    }


    if (ImGui::Button("Project Bound to closest GL(3) field")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    
    }



/*


    static ImGuiTextFilter filter;
        ImGui::Text("Filter usage:\n"
                    "  \"\"         display all lines\n"
                    "  \"xxx\"      display lines containing \"xxx\"\n"
                    "  \"xxx,yyy\"  display lines containing \"xxx\" or \"yyy\"\n"
                    "  \"-xxx\"     hide lines containing \"xxx\"");
        filter.Draw();


        
        // memset(lines, 0, sizeof(lines));
        
        // for (int i = 0; i < IM_ARRAYSIZE(lines); i++)
        //     if (filter.PassFilter(lines[i]))
        //         ImGui::BulletText("%s", lines[i]);




  if (ImGui::Button("run subroutine")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
  }

	ImGui::SameLine();
	if (ImGui::Button("hi")) {
		polyscope::warning("hi");
	}
    */



    if (ImGui::Button("Run Mint")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
        polyscope::warning("NOT IMPLEMENTED");
    }

    ImGui::SameLine();
    HelpMarker("TODO: start mint call in seperate thread");
    ImGui::SameLine();

    if (ImGui::RadioButton("Exact", cur_solver == Mint_Linear_Solver::exact))  { cur_solver = Mint_Linear_Solver::exact; } ImGui::SameLine();
    if (ImGui::RadioButton("GMRes", cur_solver == Mint_Linear_Solver::gmres))  { cur_solver = Mint_Linear_Solver::gmres; } 


 





    ImGui::ShowDemoWindow(); 

	ImGui::PopItemWidth();
}

}