#ifndef MINTGUI_H
#define MINTGUI_H

#include <Eigen/Core>
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include "polyscope/polyscope.h"

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

    enum Mint_Linear_Solver{ exact, gmres };

    enum Mint_Integrability_Mode{ free, mint };
    enum Mint_Frame_Projection{ onshell, offshell };

// TODO 
    enum Mint_Laplacian_Metric{ combinatorial, barycentric };
    enum Mint_Boundary_Conditions{ unit_normal, soft_constraint, radial, stroke };
    enum Mint_Init_Lf{ as_smooth_as_possible, random_frames, zero };



    enum Mint_Moment_Metric{ four, sec, four_plus_two, four_plus_two_tensor_two };

    enum Moments_To_Show{ second, fourth, both };
    enum Frames_To_Show{ frames, split_frames, split_moments, split_difference };



   class MintGUI
   {



    public:  
		MintGUI();

        void read_prev_state_log();
        void write_cur_state();
        void load_state_from_output_dir();


        void show_base_mesh(const std::string &id, glm::vec3 shift);
		void show_exploded_moments(Moments_To_Show moment_view_mode);
		void show_moments_2nd();
        void show_moments_4th();
        void show_moments_all();


        void show_frame_field(Frames_To_Show frame_field_view_mode);
        
        void show_gl3_frame_field(const std::string &id, glm::mat4x4 trans);
        void show_gl3_split();
        void show_integrated_quantities();

        void integrate_frame_field();
        // new file
        void convert_GL_3_frames_to_one_forms();
        void integrate_oneforms_via_spanning_tree();

        // for some reason, polyscope 
        void clear_polyscope_state();

        void rescale_structure(polyscope::Structure* m);
        void set_base_mesh();
        void set_frame_field();


        // void select_mesh();
        // void select_boundary();
        // void select_mint_output_dir();

        void set_comment();

        char* file_select_subroutine();

        

        void gui_callback();
        void gui_main_control_panel_callback();
        void gui_file_explorer_callback();   ///// reload mesh automatically.  Whole folder should use same mesh to save space
        void gui_folder_explorer_callback();
        void gui_run_mint_callback();
        void gui_file_select_BROKEN();

        void save_current_state();
        void save_polyscope_config();
        void load_prev_polyscope_config();

        FileParts fileparts(const std::string &fullpath);

        char* path_mesh;
        char* path_fra; 
        char* path_outdir;



    private:
        // matlab thread
        // json 
        char* lines;

        char* path_constraints;

        int sel_idx_mom;
        int sel_idx_fra;



        CubeCover::TetMeshConnectivity mesh;

        Eigen::MatrixXd V;
        Eigen::MatrixXi T;
        Eigen::MatrixXi F;
        Eigen::MatrixXi bdryF;



        // Moments 
        Eigen::MatrixXd M_curr;

        // Field 
        CubeCover::FrameField* field;

        // Field viz 
        Eigen::MatrixXd centroids;
        std::vector<Eigen::MatrixXd> framefieldvecs;
        Eigen::MatrixXd splitCurls;
        std::vector<Eigen::Vector2i> tree_traversal;
        std::vector<Eigen::Vector4i> tree_traversal_metadata;
        Eigen::MatrixXd treeIntegratedVals;
        float integrated_period;


        // Curve network viz stuff 
        Eigen::MatrixXd Pblack;
        Eigen::MatrixXi Eblack;
        Eigen::MatrixXd Pblue;
        Eigen::MatrixXi Eblue;
        Eigen::MatrixXd Pgreen;
        Eigen::MatrixXi Egreen;


        Mint_Linear_Solver cur_solver;
        Mint_Integrability_Mode mint_mode;
        Mint_Frame_Projection shell_mode;
        Mint_Moment_Metric metric_mode;

		Moments_To_Show moment_view_mode;
        Frames_To_Show frame_field_view_mode;
        bool showBoundary;
        bool showInteriorTets; 
        bool useSameColorRangeForAllMoments; 

        std::vector<std::string>  moment_labels; 
        std::vector<std::string>  folder_contents; 
        std::vector<std::string>  folder_contents_fra; 
        std::vector<std::string>  file_names; 
        std::vector<std::string>  file_names_fra; 
        std::vector<std::string>  file_names_mesh; 
        std::vector<std::string>  adj_folder_names; 
        std::vector<std::string>  mesh_names; // TODO

        float exploded_spacing_inp;  // fix this make user param 
        float exploded_spacing_intern;  // fix this make user param 

        float color_range_min;
        float color_range_max;

        float transparency_tets;
        float transparency_surf;

        float viz_scale;



    // bool readEdgeCurl(const std::string& edges_mintFilename, 
    //                   const std::string& edgeCurlFilename,  
    //                   Eigen::MatrixXi& edges_mint, 
    //                   Eigen::MatrixXd& edgeCurl);


    // void plot_tet_scalar_field(std::string fra_path, std::string mesh_name, std::string display_name);

    };

};

#endif