#include "ReadFancyData.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <fstream>
#include <iostream>

namespace CubeCover
{

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