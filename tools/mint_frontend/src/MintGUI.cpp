#include <igl/file_dialog_open.h>

#include "MintGUI.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

// #include <imgui/misc/cpp/imgui_stdlib.h>
// #include <polyscope/deps/imgui/imgui/imgui.h>
// #include <imgui/imgui.h>  
#include "TetMeshConnectivity.h"
#include "readMeshFixed.h"

#include <misc/cpp/imgui_stdlib.h>
#include <algorithm>
#include <cctype>
#include <string>

// This is c++14 experimental feature.  In c++17 >= this is part of STD.
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;


namespace MintFrontend
{

	MintGUI::MintGUI()
	{
        path_mesh = new char[512];
        path_constraints = new char[512];
        path_outdir = new char[512];

        
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



void MintGUI::set_base_mesh()
{
    if (!CubeCover::readMESH(path_mesh, V, T, F))
        polyscope::warning("Unable to load selected mesh");
    else{
        polyscope::registerTetMesh("tet_mesh", V, T);
        CubeCover::TetMeshConnectivity mesh(T);
        // make boundary mesh out of volume mesh
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
    
    }


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

    ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("Mesh Path", "enter rel or abs path, or use picker button below.", path_mesh, 512);
    ImGui::PopItemWidth();

    if (ImGui::Button("Pick .mesh")) {
    // executes when button is pressed
        char* tmp_path_mesh = fileSelectSubroutine();
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

    ImGui::SameLine();

    if (ImGui::Button("(re)load view mesh")) {
        // auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
        // psMesh->setTransparency(0.2);
        // psMesh->setSurfaceColor({ 0.5,0.5,0.0 });
    }

    ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("Bound Constraints (moment space)", "enter rel or abs path, or use picker button below.", path_constraints, 512);
    ImGui::PopItemWidth();



    if (ImGui::Button("Pick .bound")) {
        char* tmp_path_constraints = fileSelectSubroutine();
        strncpy(path_constraints, tmp_path_constraints, 512);
    }

    ImGui::SameLine();
    
    if (ImGui::Button("(re)Load Boundary")) {

      polyscope::warning("The chosen .bound file does not match the loaded mesh.");
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    }

     ImGui::SameLine();

    if (ImGui::Button("Project Bound to closest GL(3) field")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    }



        ImGui::PushItemWidth(300);
    ImGui::InputTextWithHint("Mint output directory", "enter rel or abs path, or use picker button below.", path_outdir, IM_ARRAYSIZE(path_outdir));
    ImGui::PopItemWidth();



    if (ImGui::Button("Pick mint output dir")) {
        char* tmp_path_outdir = fileSelectSubroutine();
        strncpy(path_outdir, tmp_path_outdir, 512);
    }

    ImGui::SameLine();
    
    if (ImGui::Button("Load dir contents")) {

      polyscope::warning("The chosen .bound file does not match the loaded mesh.");
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    }

  
    if (ImGui::Button("Project Bound to closest GL(3) field")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    }


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

    // ImGui::ShowDemoWindow(); 

	ImGui::PopItemWidth();
}

}