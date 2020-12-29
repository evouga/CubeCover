#ifndef READFRAMEFIELD_H
#define READFRAMEFIELD_H

#include <Eigen/Core>

namespace CubeCover
{
    /*
    * Reads the frame field from a .fra file, and optionally, the local
    * assignments from a .perm file. (If the permFilename is empty,
    * sets all local assignments to the identity.)
    * The |T| x 4 matrix of tetrahedral vertex indices T must be provided
    * (and readFrameField will check that the data being read from the
    * .fra and .perm files is consistent with the combinatorics of the tet
    * mesh).
    * On success, the data structures frames and assignments are populated with
    * data in the same format as expected by CubeCover() (and can be passed
    * directly in).
    *
    * Returns false if there are problems reading or parsing the input (set
    * verbose to true for diagnostic information on the console).
    */

    bool readFrameField(const std::string& fraFilename, const std::string& permFilename, const Eigen::MatrixXi& T,
        Eigen::MatrixXd& frames,
        Eigen::MatrixXi& assignments,
        bool verbose = false);

};

#endif