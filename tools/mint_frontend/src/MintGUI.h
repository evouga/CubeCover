#ifndef MINTGUI_H
#define MINTGUI_H

#include <Eigen/Core>

namespace MintFrontend
{
    /*
    * This parses the stuff dumped by dumpExtraVizStuff in the mint repo
    */

   class MintGUI
   {



    public:  
        void read_prev_state_log();
        void write_cur_state();
        void load_state_from_output_dir();

        void select_mesh();
        void select_boundary();

        void file_select_subroutine();
        void output_directory();

        void gui_callback();

        void save_current_state();


    private:
        // matlab thread
        // json 


     

    // bool readEdgeCurl(const std::string& edges_mintFilename, 
    //                   const std::string& edgeCurlFilename,  
    //                   Eigen::MatrixXi& edges_mint, 
    //                   Eigen::MatrixXd& edgeCurl);


    // void plot_tet_scalar_field(std::string fra_path, std::string mesh_name, std::string display_name);

    };

};

#endif