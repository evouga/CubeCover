#include "MIPWrapper.h"
#include <iostream>

#ifdef HAS_COMISO
#include "CoMISo/Solver/ConstrainedSolver.hh"

namespace CubeCover
{

    bool CoMISoMIPWrapper(const Eigen::SparseMatrix<double>& constraints,
        const Eigen::SparseMatrix<double>& A,
        Eigen::VectorXd& result,
        const Eigen::VectorXd& b,
        const std::vector<int>& intvars,
        double tol,
        bool verbose)
    {

        /*
     * Solves the following MIP to tolerance tol:
     * min_x  x^T A x + b x
     *  s.t.  C [x 1]^T = 0
     *        x_i is integer, for each i in intvars
     * Note that A is n x n, b is n x 1, and C is m x (n+1).
     */

        // CoMISo solvers Mx = rhs. In our case M = 2A and rhs = -b.

        // convert A to GMM matrix
        gmm::col_matrix< gmm::wsvector< double > > M(A.rows(), A.cols());
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                int row = it.row();
                int col = it.col();
                M(row, col) += 2.0 * it.value();
            }
        }

        // convert b to GMM vector
        std::vector<double> rhs(b.size());
        for (int i = 0; i < b.size(); i++)
            rhs[i] = -b[i];

        // convert constraints to GMM matrix
        gmm::row_matrix< gmm::wsvector< double > > C(constraints.rows(), constraints.cols());
        for (int k = 0; k < constraints.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(constraints, k); it; ++it) {
                int row = it.row();
                int col = it.col();
                C(row, col) += 2.0 * it.value();
            }
        }

        std::vector<double> x(A.cols());        

        COMISO::ConstrainedSolver solver;
        
        solver.misolver().set_cg_error(tol);
        solver.misolver().set_direct_rounding(false);

        // cannot be const
        std::vector<int> intvarscopy = intvars;

        solver.solve(C, M, x, rhs, intvarscopy, 0.0, verbose, false);

        result.resize(A.cols());
        for (int i = 0; i < A.cols(); i++)
            result[i] = x[i];

        return true;
    }
};

#else

namespace CubeCover
{

    bool CoMISoMIPWrapper(const Eigen::SparseMatrix<double>& constraints,
        const Eigen::SparseMatrix<double>& A,
        Eigen::VectorXd& result,
        const Eigen::VectorXd& b,
        const std::vector<int>& intvars,
        double tol,
        bool verbose)
    {
        std::cerr << "GUROBI not available" << std::endl;
        return false;
    }
};


#endif