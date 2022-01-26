#include "CurlCorrection.h"
#include "FrameField.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "TetMeshConnectivity.h"
#include <iostream>

namespace CubeCover {

    void buildCurlMatrix(int vpf, const Eigen::MatrixXd& V, const FrameField &field, Eigen::SparseMatrix<double> &C)
    {
        std::vector<Eigen::Triplet<double> > Ccoeffs;
        const TetMeshConnectivity& mesh = field.meshConnectivity();

        int nfaces = mesh.nFaces();
        int row = 0;
        for (int i = 0; i < nfaces; i++)
        {
            int t0 = mesh.faceTet(i, 0);
            int t1 = mesh.faceTet(i, 1);
            if (t0 == -1 || t1 == -1)
                continue;
            int vertids[3];
            for (int j = 0; j < 3; j++)
            {
                vertids[j] = mesh.faceVertex(i, j);
            }

            Eigen::Vector3d pts[3];
            for (int j = 0; j < 3; j++)
            {
                pts[j] = V.row(vertids[j]).transpose();
            }

            Eigen::Vector3d n = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
            n.normalize();

            Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
            AssignmentGroup o = field.faceAssignment(i);
            
            for (int j = 0; j < vpf; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {                        
                        int targetidx = o.targetVector(j);
                        double targetsign = o.targetSign(j);                        
                        Ccoeffs.push_back({ row + 3 * j + k, 3 * vpf * t0 + 3 * j + l, P(k,l) });
                        Ccoeffs.push_back({ row + 3 * j + k, 3 * vpf * t1 + 3 * targetidx + l, - targetsign * P(k,l) });
                    }
                }
            }
            row += 3 * vpf;
        }

        C.resize(row, 3 * vpf * mesh.nTets());
        C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
    }


    void curlCorrect(const Eigen::MatrixXd& V, FrameField& field, double maxCorrection, bool verbose)
    {
        int vpf = field.vectorsPerFrame();
        int ntets = field.meshConnectivity().nTets();
        Eigen::VectorXd unrolled(ntets * vpf * 3);
        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < vpf; j++)
            {
                unrolled.segment<3>(3 * vpf * i + 3 * j) = field.tetFrame(i).row(j).transpose();
            }
        }
        
        Eigen::SparseMatrix<double> C;
        buildCurlMatrix(field.vectorsPerFrame(), V, field, C);
        
        int nconstraints = C.rows();
        std::vector<Eigen::Triplet<double> > regCoeffs;
        double reg = 1e-6;
        for (int i = 0; i < nconstraints; i++)
        {
            regCoeffs.push_back({ i,i,reg });
        }
        Eigen::SparseMatrix<double> Reg(nconstraints, nconstraints);
        Reg.setFromTriplets(regCoeffs.begin(), regCoeffs.end());
        Eigen::SparseMatrix<double> Mat = Reg + C * C.transpose();
        Eigen::VectorXd rhs = C * unrolled;

        if (verbose)
        {
            std::cout << "Curl residual before: " << (C * unrolled).norm() << std::endl;
        }

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(Mat);
        Eigen::VectorXd lambda = solver.solve(rhs);

        if (verbose)
        {
            std::cout << "Curl correction solver residual: " << (Mat * lambda - rhs).norm() << std::endl;
        }

        Eigen::VectorXd delta = -C.transpose() * lambda;
        double maxmult = 0;
        for (int i = 0; i < ntets; i++)
        {
            Eigen::VectorXd frame = unrolled.segment(3 * vpf * i, 3 * vpf);
            double fnorm = frame.norm();
            Eigen::VectorXd framedelta = delta.segment(3 * vpf * i, 3 * vpf);
            double deltanorm = framedelta.norm();
            double mult = std::min(1.0, maxCorrection * fnorm / deltanorm);
            maxmult = std::max(maxmult, deltanorm / fnorm);
            unrolled.segment(3 * vpf * i, 3 * vpf) += mult * framedelta;
        }

        if (verbose)
        {
            std::cout << "Maximum correction magnitude was: " << maxmult << std::endl;
            std::cout << "Curl residual after: " << (C * unrolled).norm() << std::endl;
        }

        for (int i = 0; i < ntets; i++)
        {
            Eigen::MatrixXd newframe(vpf, 3);
            for (int j = 0; j < vpf; j++)
            {                
                newframe.row(j) = unrolled.segment<3>(3 * vpf * i + 3 * j).transpose();
            }
            field.setFrame(i, newframe);
        }
    }
};