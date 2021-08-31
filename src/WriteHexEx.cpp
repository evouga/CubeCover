#include "WriteHexEx.h"
#include <iostream>
#include <fstream>

#include <Eigen/Dense>

namespace CubeCover
{
    bool writeHexEx(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const Eigen::MatrixXd &vals)
    {
        if (4 * T.rows() != vals.rows())
        {
            std::cerr << "dimension mismatch between T and vals" << std::endl;
            return false;
        }

        std::ofstream ofs(filename);
        if (!ofs)
        {
            std::cerr << "cannot open file " << filename << " for writing" << std::endl;
            return false;
        }

        int nverts = V.rows();
        ofs << nverts << std::endl;
        for (int i = 0; i < nverts; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (j != 0)
                    ofs << " ";
                ofs << V(i, j);
            }
            ofs << std::endl;
        }

        int ntets = T.rows();
        int nvals = vals.cols();
        if (nvals != 3)
        {
            std::cerr << "warning: number of parameter values per vertex is not equal to 3. HexEx file format only officially supports 3D parameterization" << std::endl;
        }
        ofs << ntets << std::endl;

        for (int i = 0; i < ntets; i++)
        {
            bool hasVol = true; 
            // for (int dir = 0; dir < 3; dir++)
            // {            
            //     int sample1 = (dir + 1) % 3;
            //     int sample2 = (dir + 2) % 3;

            //     for (int vert1 = 0; vert1 < 4; vert1++)
            //     {
            //         for (int vert2 = vert1 + 1; vert2 < 4; vert2++)
            //         {
            //             for (int vert3 = vert2 + 1; vert3 < 4; vert3++)
            //             {
            //                 double s11val = vals(4 * i + vert1, sample1);
            //                 double s12val = vals(4 * i + vert2, sample1);
            //                 double s13val = vals(4 * i + vert3, sample1);
            //                 double s21val = vals(4 * i + vert1, sample2);
            //                 double s22val = vals(4 * i + vert2, sample2);
            //                 double s23val = vals(4 * i + vert3, sample2);
                                    
            //                 Eigen::Matrix2d M;
            //                 M << s12val - s11val, s13val - s11val,
            //                      s22val - s21val, s23val - s21val;

            //                 if (M.determinant() == 0)
            //                 {
            //                     hasVol = false;
            //                 }
            //             }
            //         }
            //     }
            // }

            if(hasVol)
            {
                // int tmp = T(i, 0);
                // T(i, 0) = T(i, 1);
                // T(i, 1) = tmp;
                // Eigen::VectorXd tmpv = vals.row(4*i);
                // vals.row(4*i) = vals.row(4*i + 1);
                // vals.row(4*i + 1) = tmpv;
                for (int j = 0; j < 4; j++)
                {
                    if (j != 0)
                        ofs << " ";
                    ofs << T(i, j);
                }
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < nvals; k++)
                    {
                        ofs << " " << vals(4 * i + j, k);
                    }
                }
                ofs << std::endl;
            }


        }
        if (!ofs)
        {
            std::cerr << "error writing to file " << filename << std::endl;
            return false;
        }
        return true;
    }
};