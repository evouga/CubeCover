#ifndef MIPWRAPPER_H
#define MIPWRAPPER_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include "CubeCover.h"

namespace CubeCover {

    /*
     * Solves the following MIP to tolerance tol:
     * min_x  x^T A x + b x
     *  s.t.  C [x 1]^T = 0
     *        x_i is integer, for each i in intvars
     * Note that A is n x n, b is n x 1, and C is m x (n+1).
     */

    bool GurobiMIPWrapper(const Eigen::SparseMatrix<double>& constraints,
        const Eigen::SparseMatrix<double>& A,
        Eigen::VectorXd& result,
        const Eigen::VectorXd& b,
        const std::vector<int>& intvars,
        double tol,
        bool verbose);

    bool CoMISoMIPWrapper(const Eigen::SparseMatrix<double>& constraints,
        const Eigen::SparseMatrix<double>& A,
        Eigen::VectorXd& result,
        const Eigen::VectorXd& b,
        const std::vector<int>& intvars,
        double tol,
        bool verbose);


    inline bool MIP(const Eigen::SparseMatrix<double>& constraints,
        const Eigen::SparseMatrix<double>& A,
        Eigen::VectorXd& result,
        const Eigen::VectorXd& b,
        const std::vector<int>& intvars,
        const CubeCover::CubeCoverOptions& opts)
    {
        if (opts.solver == CubeCoverOptions::MIPSolver::MS_GUROBI)
        {
            return GurobiMIPWrapper(constraints, A, result, b, intvars, opts.MIPtol, opts.verbose);
        }
        else if (opts.solver == CubeCoverOptions::MIPSolver::MS_COMISO)
        {
            return CoMISoMIPWrapper(constraints, A, result, b, intvars, opts.MIPtol, opts.verbose);
        }

        // unreachable
        return false;
    }
    
};

#endif
