#include <igl/file_dialog_open.h>

#include "MintGUI.h"
#include "polyscope/polyscope.h"

// #include <imgui/misc/cpp/imgui_stdlib.h>
// #include <polyscope/deps/imgui/imgui/imgui.h>
// #include <imgui/imgui.h>  

#include <misc/cpp/imgui_stdlib.h>

namespace MintFrontend
{

	MintGUI::MintGUI()
	{

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

    char path_mesh[512] = "";
    char path_constraints[512] = "";
    char path_outdir[512] = "";

    // std::cout << rightWindowsWidth << std::endl;

    ImGui::PushItemWidth(512);
    ImGui::InputTextWithHint("Mesh Path", "enter rel or abs path, or use picker button below.", path_mesh, IM_ARRAYSIZE(path_mesh));
    ImGui::PopItemWidth();

    if (ImGui::Button("Pick .mesh")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    }

    ImGui::SameLine();

    if (ImGui::Button("(re)load view mesh")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
    }

    ImGui::PushItemWidth(512);
    ImGui::InputTextWithHint("Bound Constraints (moment space)", "enter rel or abs path, or use picker button below.", path_constraints, IM_ARRAYSIZE(path_constraints));
    ImGui::PopItemWidth();



    if (ImGui::Button("Pick .bound")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
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



        ImGui::PushItemWidth(512);
    ImGui::InputTextWithHint("Mint output directory", "enter rel or abs path, or use picker button below.", path_outdir, IM_ARRAYSIZE(path_outdir));
    ImGui::PopItemWidth();



    if (ImGui::Button("Pick mint output dir")) {
    // executes when button is pressed
    // directory_path = fileSelectSubroutine();
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


        
        memset(lines, 0, sizeof(lines));
        
        for (int i = 0; i < IM_ARRAYSIZE(lines); i++)
            if (filter.PassFilter(lines[i]))
                ImGui::BulletText("%s", lines[i]);




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