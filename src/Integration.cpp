#include "Integration.h"
#include "GurobiMIPWrapper.h"
#include <iostream>
#include "CutToSimplyConnected.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <set>
#include <Eigen/Dense>
#include <vector>

namespace CubeCover {

    class UnionFind {
    public:
        UnionFind(int nelements)
        {
            parent.resize(nelements);
            sign.resize(nelements);
            for (int i = 0; i < nelements; i++)
            {
                parent[i] = i;
                sign[i] = 1;
            }
        }

        void find(int idx, int &label, int &labelsign)
        {
            int root = idx;
            int signtoroot = 1;
            while (parent[root] != root)
            {
                signtoroot *= sign[root];
                root = parent[root];
            }
            labelsign = signtoroot;          
            while (parent[idx] != root)
            {
                int p = parent[idx];
                int s = signtoroot * sign[idx];
                parent[idx] = root;
                sign[idx] = signtoroot;
                idx = p;
                signtoroot = s;
            }
            label = root;              
        }

        void unionTogether(int idx1, int idx2, int labelsign)
        {
            int l1, ls1, l2, ls2;
            find(idx1, l1, ls1);
            find(idx2, l2, ls2);
            parent[l1] = l2;
            sign[l1] = labelsign * ls1 * ls2;
        }

    private:
        std::vector<int> parent;
        std::vector<int> sign;
    };


    bool integrate(const Eigen::MatrixXd& V, const FrameField& field, Eigen::MatrixXd& soupValues, double scale, 
        double MIPtol, CubeCoverOptions::ParameterizationType type, bool forceBoundaryAlignment, bool verbose)
    {
        if (type == CubeCoverOptions::ParameterizationType::PT_INTEGERGRID)
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

        soupValues.resize(4 * ntets, vpf);
        if (ntets == 0)
            return true;

        std::set<int> cutfaceset;
        for (auto it : cutFaces)
            cutfaceset.insert(it);

        // three parameter values per soup vertex
        int soupdofs = 4 * ntets * vpf;
        UnionFind uf(soupdofs);

        // merge together equivalent DOFs
        int nfaces = mesh.nFaces();
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
                    if (!cutfaceset.count(i))
                    {
                        uf.unionTogether(fromidx[j] + k, toidx[j] + destindex, destsign);
                    }                    
                }
            }            
        }

        std::map<int, int> labelmap;
        int reduceddofs = 0;
        for (int i = 0; i < soupdofs; i++)
        {
            int l, s;
            uf.find(i, l, s);
            auto it = labelmap.find(l);
            if (it == labelmap.end())
            {
                labelmap[l] = reduceddofs++;
            }
        }

        // three jump variables across each face
        int jumpdofs = vpf * cutFaces.size();



        if(verbose)
            std::cout << "MIP problem has " << reduceddofs + jumpdofs << " total variables and " << jumpdofs << " integer variables" << std::endl;

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
                    if (cutfaceset.count(i))
                    {
                        int froml, froms;
                        uf.find(fromidx[j] + k, froml, froms);
                        Ccoeffs.push_back({ nconstraints, labelmap[froml], double(froms) });
                        int tol, tos;
                        uf.find(toidx[j] + destindex, tol, tos);                        
                        Ccoeffs.push_back({ nconstraints, labelmap[tol], double(-tos * destsign) });
                        Ccoeffs.push_back({ nconstraints, reduceddofs + intdofidx + k, 1.0 });
                        nconstraints++;
                    }
                    
                }
            }
            if (cutfaceset.count(i))
                intdofidx += vpf;
        }
        // pin one corner
        for (int i = 0; i < vpf; i++)
        {
            int l, s;
            uf.find(i, l, s);
            Ccoeffs.push_back({ nconstraints, labelmap[l], 1.0 });
            nconstraints++;
        }
        if(verbose)
            std::cout << "Enforcing " << nconstraints << " constraints" << std::endl;
        assert(intdofidx == jumpdofs);
        Eigen::SparseMatrix<double> C(nconstraints, reduceddofs + jumpdofs + 1);
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
                        int label, sign;
                        uf.find(12 * i + vpf * j + l, label, sign);
                        Dcoeffs.push_back({ 3 * vpf * i + 3 * l + k, labelmap[label], sign * Grad(k,j) });
                    }
                }
            }
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < vpf; l++)
                {
                    rhs[3 * vpf * i + 3 * l + k] = scale * field.tetFrame(i)(l, k);
                }
            }
        }
        Eigen::SparseMatrix<double> D(3 * vpf * ntets, reduceddofs + jumpdofs);
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
        std::vector<int> intdofs;
        
        if (type == CubeCoverOptions::ParameterizationType::PT_SEAMLESS)
        {
            intdofs.resize(jumpdofs);
            for (int i = 0; i < jumpdofs; i++)
            {
                intdofs[i] = reduceddofs + i;
            }
        }

        if (!GurobiMIPWrapper(C, MIPA, result, MIPrhs, intdofs, 1e-6, verbose))
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
                    std::cout << result[reduceddofs + i] << " ";
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
                        int label, sign;
                        uf.find(vpf * 4 * i + vpf * j + k, label, sign);
                        soupValues(4 * i + j, k) = result[labelmap[label]] * sign;
                    }
                }
            }
        }
        return true;
    }
    
};