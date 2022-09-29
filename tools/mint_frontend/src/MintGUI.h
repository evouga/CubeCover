#ifndef MINTGUI_H
#define MINTGUI_H

#include <Eigen/Core>
#include "TetMeshConnectivity.h"

namespace MintFrontend
{
    /*
    * This parses the stuff dumped by dumpExtraVizStuff in the mint repo
    */

    struct FileParts
    {
        std::string path; //!< containing folder, if provided, including trailing slash
        std::string name; //!< base file name, without extension
        std::string ext;  //!< extension, including '.'
    };


   class MintGUI
   {



    public:  
        void read_prev_state_log();
        void write_cur_state();
        void load_state_from_output_dir();


        void set_base_mesh();
        // void select_mesh();
        // void select_boundary();
        // void select_mint_output_dir();

        void set_comment();

        char* file_select_subroutine();
        

        void gui_callback();

        void save_current_state();

        FileParts fileparts(const std::string &fullpath);
        

		MintGUI();

    private:
        // matlab thread
        // json 
        char* lines;
        char* path_mesh;
        char* path_constraints;
        char* path_outdir;


        CubeCover::TetMeshConnectivity mesh;

        Eigen::MatrixXd V;
        Eigen::MatrixXi T;
        Eigen::MatrixXi F;
        Eigen::MatrixXi bdryF;


     

    // bool readEdgeCurl(const std::string& edges_mintFilename, 
    //                   const std::string& edgeCurlFilename,  
    //                   Eigen::MatrixXi& edges_mint, 
    //                   Eigen::MatrixXd& edgeCurl);


    // void plot_tet_scalar_field(std::string fra_path, std::string mesh_name, std::string display_name);

    };

};

#endif