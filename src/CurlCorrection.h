#ifndef CURLCORRECTION_H
#define CURLCORRECTION_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class FrameField;

namespace CubeCover
{

	void buildCurlMatrix(int vpf, const Eigen::MatrixXd& V, const FrameField &field, Eigen::SparseMatrix<double> &C);
    
    void curlCorrect(const Eigen::MatrixXd& V, FrameField& field, double maxCorrection, bool verbose);
};

#endif