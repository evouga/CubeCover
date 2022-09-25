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
        void select_mint_output_dir();

        void set_comment();

        void file_select_subroutine();
        

        void gui_callback();

        void save_current_state();

		MintGUI();

    private:
        // matlab thread
        // json 
        char* lines[];


     

    // bool readEdgeCurl(const std::string& edges_mintFilename, 
    //                   const std::string& edgeCurlFilename,  
    //                   Eigen::MatrixXi& edges_mint, 
    //                   Eigen::MatrixXd& edgeCurl);


    // void plot_tet_scalar_field(std::string fra_path, std::string mesh_name, std::string display_name);

    };

};

#endif