#include "ReadHexEx.h"
#include <iostream>
#include <fstream>

namespace CubeCover
{
    bool readHexEx(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& T, Eigen::MatrixXd &vals)
    {
        std::ifstream ifs(filename);
        if (!ifs)
        {
            std::cerr << "cannot open file " << filename << " for reading" << std::endl;
            return false;
        }

        int nverts;
        ifs >> nverts;
        if (!ifs)
            return false;

        V.resize(nverts, 3);
        for (int i = 0; i < nverts; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                ifs >> V(i, j);
            }            
        }

        int ntets;
        ifs >> ntets;
        if (!ifs)
            return false;

        T.resize(ntets, 4);
        vals.resize(4 * ntets, 3);

        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                ifs >> T(i, j);
            }
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ifs >> vals(4 * i + j, k);
                }
            }            
        }
        if (!ifs)
        {
            return false;
        }
        return true;
    }
};