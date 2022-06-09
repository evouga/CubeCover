#include "ReadFancyData.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <fstream>
#include <iostream>
#include <Eigen/Core>

#include "polyscope/polyscope.h"
#include "polyscope/volume_mesh.h"

namespace CubeCover
{

    void plot_tet_scalar_field(std::string fra_path, std::string mesh_name, std::string display_name)
    {
        std::ifstream ifs(fra_path);
        Eigen::VectorXd tet_values;
        int ntets;
        ifs >> ntets;
        tet_values.resize(ntets);
        
        for (int i = 0; i < ntets; i++)
        {
            ifs >> tet_values(i);
        }
        ifs.close();
        Eigen::VectorXd tmp = tet_values;
        std::sort(tmp.data(), tmp.data() + tmp.size() ) ;
        std::pair<double,double> init_colormap_range = std::pair<double,double>(tmp(0), tmp(int( tmp.rows() * .9 ) ));

        auto colorQ = polyscope::getVolumeMesh("tet_mesh")->addCellScalarQuantity(display_name, tet_values);
        colorQ->setEnabled(false);
        colorQ->setMapRange(init_colormap_range);

        return;


    }


    bool readEdgeCurl(const std::string& edges_mintFilename, const std::string& edgeCurlFilename, Eigen::MatrixXi& edgedata_mint, Eigen::MatrixXd& edgeCurl)
    {
        std::ifstream ifs(edgeCurlFilename);
        std::ifstream ifs_edgedata(edges_mintFilename);
        
        bool verbose = true;
        if (!ifs)
        {
            if (verbose)
                std::cerr << "Cannot open file: " << edgeCurlFilename << std::endl;
            return false;
        }

        int nedges;
        ifs >> nedges;

        edgeCurl.resize(nedges, 1);
        for ( int i = 0; i < nedges; i++)
        {
            ifs >> edgeCurl(i);
            if (!ifs)
            {
                if (verbose)
                    std::cerr << "Error reading curl data" << std::endl;
                return false;
            
            }
        }

        edgedata_mint.resize(nedges, 4);
        for ( int i = 0; i < nedges; i++)
        {
            for (int j = 0; j < 4; j ++)
            {
                ifs_edgedata >> edgedata_mint(i, j);
                if (!ifs_edgedata)
                {
                    if (verbose)
                        std::cerr << "Error reading edge data" << std::endl;
                    return false;
                
                }
            }
        }


        ifs.close();
        ifs_edgedata.close();

        return true;

    }

};