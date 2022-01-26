#ifndef READFANCYDATA_H
#define READFANCYDATA_H

#include <Eigen/Core>

namespace CubeCover
{
    /*
    * This parses the stuff dumped by dumpExtraVizStuff in the mint repo
    */

    bool readEdgeCurl(const std::string& edges_mintFilename, 
                      const std::string& edgeCurlFilename,  
                      Eigen::MatrixXi& edges_mint, 
                      Eigen::MatrixXd& edgeCurl);

};

#endif