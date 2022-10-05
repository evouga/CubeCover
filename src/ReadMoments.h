#ifndef READMOMENTS_H
#define READMOMENTS_H

#include <string>
#include <Eigen/Core>

namespace CubeCover
{
    /*
     * Reads (currently) 4th and 2nd moments of a tetrahedral mesh and boundary .mom file. Only 
     * supports 3D parameterizations. Follows the file format specification
     * of the mint library.
     * Returns true if loading is successful. On success, ME will be a
     * matrix of (ntets+nboundelements) x 22 of moments per interior volume element and surface element
     */
    bool readMoments(const std::string& filename, Eigen::MatrixXd& ME, bool verbose = false);
};

#endif