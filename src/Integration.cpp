#include "Integration.h"
#include "GurobiMIPWrapper.h"
#include <iostream>
#include "CutToSimplyConnected.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <set>
#include <Eigen/Dense>

namespace CubeCover {

    bool integrate(const Eigen::MatrixXd& V, const FrameField& field, Eigen::MatrixXd& soupValues, double scale, 
        double MIPtol, bool integerGrid, bool forceBoundaryAlignment, bool verbose)
    {
        if (integerGrid)
        {
            std::cerr << "Integer-grid parameterization not yet supported" << std::endl;
            return false;
        }
        if (forceBoundaryAlignment)
        {
            std::cerr << "Boundary alignment not yet supported" << std::endl;
            return false;
        }

        if (verbose)
        {
            std::cout << "Cutting to a simply-connected domain..." << std::endl;
        }

        const TetMeshConnectivity& mesh = field.meshConnectivity();
        int nsingular = field.nSingularEdges();
        std::vector<int> singularEdges;
        for (int i = 0; i < nsingular; i++)
            singularEdges.push_back(field.singularEdge(i));

        std::vector<int> cutFaces;
        cutToSimplyConnected(mesh, singularEdges, cutFaces, true);
        if (verbose)
        {
            std::cout << "...done, found " << cutFaces.size() << " seam faces" << std::endl;
        }

        int ntets = mesh.nTets();

        int vpf = field.vectorsPerFrame();

        // three parameter values per soup vertex
        int regulardofs = vpf * 4 * ntets;

        soupValues.resize(4 * ntets, vpf);
        if (ntets == 0)
            return true;

        // three jump variables across each face
        int jumpdofs = vpf * cutFaces.size();

        std::set<int> cutfaceset;
        for (auto it : cutFaces)
            cutfaceset.insert(it);

        if(verbose)
            std::cout << "MIP problem has " << regulardofs << " real variables and " << jumpdofs << " integer variables" << std::endl;

        int nconstraints = 0;
        int intdofidx = 0;
        std::vector<Eigen::Triplet<double> > Ccoeffs;
        for (int i = 0; i < mesh.nFaces(); i++)
        {
            if (mesh.isBoundaryFace(i))
                continue;

            int fromtet = mesh.faceTet(i, 0);
            int totet = mesh.faceTet(i, 1);
            AssignmentGroup o = field.faceAssignment(i);
            int fromidx[3];
            int toidx[3];
            for (int j = 0; j < 3; j++)
            {
                fromidx[j] = vpf * 4 * fromtet + vpf * mesh.faceTetVertexIndex(i, 0, j);
                toidx[j] = vpf * 4 * totet + vpf * mesh.faceTetVertexIndex(i, 1, j);
            }
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < vpf; k++)
                {
                    int destindex = o.targetVector(k);
                    int destsign = o.targetSign(k);
                    Ccoeffs.push_back({ nconstraints, fromidx[j] + k, 1.0 });
                    Ccoeffs.push_back({ nconstraints, toidx[j] + destindex, -double(destsign) });
                    if (cutfaceset.count(i))
                    {
                        Ccoeffs.push_back({ nconstraints, regulardofs + intdofidx + k, 1.0 });
                    }
                    nconstraints++;
                }
            }
            if (cutfaceset.count(i))
                intdofidx += vpf;
        }
        // pin one corner
        for (int i = 0; i < vpf; i++)
        {
            Ccoeffs.push_back({ nconstraints, i, 1.0 });
            nconstraints++;
        }
        if(verbose)
            std::cout << "Enforcing " << nconstraints << " constraints" << std::endl;
        assert(intdofidx == jumpdofs);
        Eigen::SparseMatrix<double> C(nconstraints, regulardofs + jumpdofs + 1);
        C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());



        // objective function
        // D: computes a per-tet gradient of each phi function
        std::vector<Eigen::Triplet<double> > Dcoeffs;
        Eigen::VectorXd rhs(3 * vpf * ntets);
        for (int i = 0; i < ntets; i++)
        {
            Eigen::Matrix3d B;
            for (int j = 0; j < 3; j++)
            {
                B.row(j) = V.row(mesh.tetVertex(i, j + 1)) - V.row(mesh.tetVertex(i, 0));
            }
            Eigen::Matrix<double, 3, 4> S;
            S << -1, 1, 0, 0,
                -1, 0, 1, 0,
                -1, 0, 0, 1;
            Eigen::Matrix<double, 3, 4> Grad = B.inverse() * S;
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < vpf; l++)
                    {
                        Dcoeffs.push_back({ 3 * vpf * i + 3 * l + k, 12 * i + vpf * j + l, Grad(k,j) });
                    }
                }
            }
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < vpf; l++)
                {
                    rhs[3 * vpf * i + 3 * l + k] = scale * field.faceFrame(i)(l, k);
                }
            }
        }
        Eigen::SparseMatrix<double> D(3 * vpf * ntets, regulardofs + jumpdofs);
        D.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());

        // inner product
        std::vector<Eigen::Triplet<double> > Mcoeffs;
        double totvol = 0;
        for (int i = 0; i < ntets; i++)
        {
            Eigen::Vector3d e1 = V.row(mesh.tetVertex(i, 1)) - V.row(mesh.tetVertex(i, 0));
            Eigen::Vector3d e2 = V.row(mesh.tetVertex(i, 2)) - V.row(mesh.tetVertex(i, 0));
            Eigen::Vector3d e3 = V.row(mesh.tetVertex(i, 3)) - V.row(mesh.tetVertex(i, 0));
            double vol = 1.0 / 6.0 * std::fabs(e1.cross(e2).dot(e3));
            totvol += vol;
            for (int j = 0; j < 3 * vpf; j++)
            {
                Mcoeffs.push_back({ 3 * vpf * i + j, 3 * vpf * i + j, vol });
            }
        }
        Eigen::SparseMatrix<double> M(3 * vpf * ntets, 3 * vpf * ntets);
        M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());

        M /= totvol;

        Eigen::SparseMatrix<double> fitterm = 0.5 * D.transpose() * M * D; // 12T x 12T

        Eigen::SparseMatrix<double> MIPA = fitterm;
        Eigen::VectorXd MIPrhs = -D.transpose() * M * rhs;

        Eigen::VectorXd result;
        std::vector<int> intdofs(jumpdofs);
        for (int i = 0; i < jumpdofs; i++)
        {
            intdofs[i] = regulardofs + i;
        }        

        if (!GurobiMIPWrapper(C, MIPA, result, MIPrhs, intdofs, 1e-6))
        {
            if (verbose)
            {
                std::cerr << "Problem with MIP solve!" << std::endl;
            }
            soupValues.resize(4 * ntets, vpf);
            soupValues.setZero();
        }
        else
        {
            if (verbose)
            {
                std::cout << "Integer variables: ";
                for (int i = 0; i < jumpdofs; i++)
                    std::cout << result[regulardofs + i] << " ";
                std::cout << std::endl;
            }
            if (verbose)
            {
                std::cout << "Problem residual: " << (MIPA * result - MIPrhs).norm() << std::endl;
            }

            soupValues.resize(4 * ntets, vpf);
            for (int i = 0; i < ntets; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < vpf; k++)
                    {
                        soupValues(4 * i + j, k) = result[vpf * 4 * i + vpf * j + k];
                    }
                }
            }
        }
        return true;
    }
    
};