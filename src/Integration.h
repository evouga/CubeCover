#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <Eigen/Core>
#include "CubeCover.h"

namespace CubeCover {

    class FrameField;

    /*
     * Performs the actual integration, by solving the problem
     *  min_phi \sum_{tets i} \sum_{frame vector j} ||\nabla \phi_j - f^i_j||^2
     *      s.t. \phi differing by an integer jump, as well as the provided 
     *           local reassignment of vectors and signs, across each interior
     *           face.
     * using the Gurobi MIP solver.
     * Optionally also enforces integer-grid and boundary alignment constraints.
     */
    bool integrate(const Eigen::MatrixXd& V, const FrameField& field, Eigen::MatrixXd& soupValues, const CubeCoverOptions &options);

};

#endif