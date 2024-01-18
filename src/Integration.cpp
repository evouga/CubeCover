#include "Integration.h"
#include "MIPWrapper.h"
#include <iostream>
#include "CutToSimplyConnected.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <set>
#include <Eigen/Dense>
#include <vector>

namespace CubeCover {

    static Eigen::Vector3d faceNormal(const Eigen::MatrixXd& V, const TetMeshConnectivity& mesh, int face)
    {
        Eigen::Vector3d p0 = V.row(mesh.faceVertex(face, 0)).transpose();
        Eigen::Vector3d p1 = V.row(mesh.faceVertex(face, 1)).transpose();
        Eigen::Vector3d p2 = V.row(mesh.faceVertex(face, 2)).transpose();
        Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
        n /= n.norm();
        return n;
    }

    static int closestFrameVec(const Eigen::Vector3d& dir, const Eigen::MatrixXd& frame)
    {
        int nvecs = frame.rows();
        int bestrow = -1;
        double bestdot = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < nvecs; i++)
        {
            Eigen::Vector3d framevec = frame.row(i).transpose();
            double dot = std::fabs(dir.dot(framevec) / framevec.norm());
            if (dot > bestdot)
            {
                bestrow = i;
                bestdot = dot;
            }
        }
        return bestrow;
    }

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


    bool integrate(const Eigen::MatrixXd& V, const FrameField& field, Eigen::MatrixXd& soupValues, const CubeCoverOptions &opt)
    {
        bool verbose = opt.verbose;

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
        
        // merge together equivalent DOFs
        int nfaces = mesh.nFaces();
        UnionFind uf(soupdofs);

        for (int i = 0; i < mesh.nFaces(); i++)
        {
            if (mesh.isBoundaryFace(i) || cutfaceset.count(i))
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
                    uf.unionTogether(fromidx[j] + k, toidx[j] + destindex, destsign);                    
                }
            }            
        }

        // boundary faces must lie on the same integer plane if alignment is on
        if (opt.boundaryConditions == CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER)
        {
            for (int i = 0; i < mesh.nFaces(); i++)
            {
                if (mesh.isBoundaryFace(i))
                {
                    int tet = mesh.faceTet(i, 0);
                    if (tet == -1)
                        tet = mesh.faceTet(i, 1);
                    assert(tet != -1);

                    Eigen::MatrixXd bdryframe = field.tetFrame(tet);
                    Eigen::Vector3d normal = faceNormal(V, mesh, i);
                    int alignedvec = closestFrameVec(normal, bdryframe);
                    assert(alignedvec != -1);
                    std::set<int> faceverts;
                    for (int j = 0; j < 3; j++)
                    {
                        faceverts.insert(mesh.faceVertex(i, j));
                    }
                    std::vector<int> tofuse;
                    for (int j = 0; j < 4; j++)
                    {
                        if (faceverts.count(mesh.tetVertex(tet, j)))
                        {
                            int soupidx = 4 * vpf * tet + vpf * j + alignedvec;
                            tofuse.push_back(soupidx);
                        }
                    }
                    for (int j = 1; j < tofuse.size(); j++)
                    {
                        uf.unionTogether(tofuse[0], tofuse[j], 1);
                    }
                }
            }
        }

        std::set<int> singularcurvesoupdofs;

        if (opt.parameterizationType == CubeCoverOptions::ParameterizationType::PT_INTEGERGRID)
        {
            int nsingedges = field.nSingularEdges();
            for (int i = 0; i < nsingedges; i++)
            {
                int edge = field.singularEdge(i);
                assert(!mesh.isBoundaryEdge(edge));
                int nbtets = mesh.nEdgeTets(edge);
                int vs[2];
                vs[0] = mesh.edgeVertex(edge, 0);
                vs[1] = mesh.edgeVertex(edge, 1);

                // must check all neighboring tets because there might be multiple connected components of tets separated by seam faces

                for (int j = 0; j < nbtets; j++)
                {
                    // circulate, starting from here
                    AssignmentGroup o(vpf);
                    for (int k = 0; k < nbtets; k++)
                    {
                        int nb = (j + k) % nbtets;
                        int tet = mesh.edgeTet(edge, nb);
                        int faceidx = mesh.edgeTetFaceIndex(edge, nb, 1);
                        int face = mesh.tetFace(tet, faceidx);
                        int orient = mesh.tetFaceOrientation(tet, faceidx);
                        AssignmentGroup faceo = field.faceAssignment(face);

                        if (orient == 1)
                            faceo = faceo.inverse();
                        o = faceo * o;
                    }
                    for (int k = 0; k < vpf; k++)
                    {
                        if (o.targetVector(k) != k || o.targetSign(k) != 1)
                        {
                            int dofs[2];
                            for (int l = 0; l < 2; l++)
                            {
                                int vertidx = -1;
                                int tet = mesh.edgeTet(edge, j);
                                for (int m = 0; m < 4; m++)
                                {
                                    if (mesh.tetVertex(tet, m) == vs[l])
                                        vertidx = m;
                                }
                                assert(vertidx != -1);
                                dofs[l] = 4 * tet * vpf + vpf * vertidx + k;                                
                            }
                            uf.unionTogether(dofs[0], dofs[1], 1);
                            singularcurvesoupdofs.insert(dofs[0]);
                            singularcurvesoupdofs.insert(dofs[1]);
                        }
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

        struct Constraint
        {
            Constraint(int dof1, int sign1,
                int dof2, int sign2,
                int intdof, int sign3)
            {
                // sort into canonical form
                // dof1 < dof2 and dof1 positive
                if (dof1 > dof2)
                {
                    std::swap(dof1, dof2);
                    std::swap(sign1, sign2);
                }
                if (sign1 == -1)
                {
                    sign2 *= -1;
                    sign3 *= -1;
                }
                d1 = dof1;
                d2 = dof2;
                intd = intdof;
                sd2 = sign2;
                sintd = sign3;
            }
            int d1, d2, intd;
            int sd2, sintd;
        };


        std::vector<Constraint> constraints;
        for (int i = 0; i < mesh.nFaces(); i++)
        {                        
            if(cutfaceset.count(i))
            {
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
                        int froml, froms;
                        uf.find(fromidx[j] + k, froml, froms);                        
                        int tol, tos;
                        uf.find(toidx[j] + destindex, tol, tos);                        
                        constraints.push_back({ labelmap[froml], froms, labelmap[tol], -tos * destsign, vpf * i + k, 1 });                        
                    }
                }
            }
        }

        UnionFind cuf(nfaces* vpf);
        std::set<int> implicitzeros;

        // group by first vertex
        int totalconstraints = constraints.size();
        std::map<int, std::vector<int> > constraintsbyvert;
        std::map<int, std::map<std::pair<int, int>, int> > representatives;
        for (int i = 0; i < totalconstraints; i++)
        {
            constraintsbyvert[constraints[i].d1].push_back(i);
        }
        for (auto& it : constraintsbyvert)
        {
            std::map<std::pair<int, int>, int> &representativeconstraints = representatives[it.first];
            for (auto& it2 : it.second)
            {
                int d2 = constraints[it2].d2;
                int s2 = constraints[it2].sd2;
                auto prev = representativeconstraints.find({ d2,s2 });
                if (prev == representativeconstraints.end())
                {
                    representativeconstraints[{d2, s2}] = it2;
                }
                else
                {
                    int previntd = constraints[prev->second].intd;
                    int prevsintd = constraints[prev->second].sintd;
                    int intd = constraints[it2].intd;
                    int sintd = constraints[it2].sintd;

                    int prevlabel, prevsign;
                    int label, sign;
                    cuf.find(previntd, prevlabel, prevsign);
                    cuf.find(intd, label, sign);
                    if (prevlabel == label && sign*sintd != prevsign*prevsintd)
                    {
                        implicitzeros.insert(intd);
                    }
                    cuf.unionTogether(intd, previntd, sintd * prevsintd);                    
                }
            }
        }

        std::set<int> implicitzerolabels;
        for (auto it : implicitzeros)
        {
            int label, sign;
            cuf.find(it, label, sign);
            implicitzerolabels.insert(label);
        }
        int jumpdofs = 0;
        std::map<int, int> jumplabelmap;        
        for (auto& it : representatives)
        {
            for (auto& it2 : it.second)
            {
                int itdof = constraints[it2.second].intd;
                int label, sign;
                cuf.find(itdof, label, sign);
                if (!implicitzerolabels.count(label))
                {
                    auto idx = jumplabelmap.find(label);
                    if (idx == jumplabelmap.end())
                    {
                        jumplabelmap[label] = jumpdofs++;
                    }
                }
            }
        }

        if (verbose)
        {
            std::cout << "MIP problem has " << reduceddofs << " parameter variables and " << jumpdofs << " transition jumps" << std::endl;
            std::cout << "(" << implicitzerolabels.size() << " transition jumps implicitly zero)" << std::endl;
        }

        std::set<int> integerreduceddofs;

        
        int nconstraints = 0;
        std::vector<Eigen::Triplet<double> > Ccoeffs;
        for (auto& it : representatives)
        {
            for (auto& it2 : it.second)
            {
                int intdof = constraints[it2.second].intd;
                int intdsign = constraints[it2.second].sintd;

                int froml = it.first;
                int tol = it2.first.first;
                int tos = it2.first.second;
                
                Ccoeffs.push_back({ nconstraints, froml, 1.0 });
                Ccoeffs.push_back({ nconstraints, tol, double(tos) });
                int idlabel, idsign;
                cuf.find(intdof, idlabel, idsign);
                if (!implicitzerolabels.count(idlabel))
                {
                    Ccoeffs.push_back({ nconstraints, reduceddofs + jumplabelmap[idlabel], double(idsign*intdsign) });
                }
                nconstraints++;
            }
        }

        bool forceBoundaryAlignment = (opt.boundaryConditions == CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER);

        for (int i = 0; i < mesh.nFaces(); i++)
        {            
            if (mesh.isBoundaryFace(i) && forceBoundaryAlignment)
            {                
                int tet = mesh.faceTet(i, 0);
                if (tet == -1)
                    tet = mesh.faceTet(i, 1);
                assert(tet != -1);

                Eigen::MatrixXd bdryframe = field.tetFrame(tet);
                Eigen::Vector3d normal = faceNormal(V, mesh, i);
                int alignedvec = closestFrameVec(normal, bdryframe);
                assert(alignedvec != -1);
                std::set<int> faceverts;
                for (int j = 0; j < 3; j++)
                {
                    faceverts.insert(mesh.faceVertex(i, j));
                }
                for (int j = 0; j < 4; j++)
                {
                    if (faceverts.count(mesh.tetVertex(tet, j)))
                    {
                        int soupidx = 4 * vpf * tet + vpf * j + alignedvec;
                        int l, s;
                        uf.find(soupidx, l, s);
                        integerreduceddofs.insert(labelmap[l]);
                    }
                }
            }            
        }

        if (opt.parameterizationType == CubeCoverOptions::ParameterizationType::PT_INTEGERGRID)
        {
            for (auto it : singularcurvesoupdofs)
            {
                int label, sign;
                uf.find(it, label, sign);
                integerreduceddofs.insert(labelmap[label]);
            }
        }


        if (verbose)
        {
            std::cout << integerreduceddofs.size() << " parameter DOFs are integer" << std::endl;
        }

        // pin one corner
        int maxintcoords = -1;
        int pindof = 0;
        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int intcoords = 0;
                for (int k = 0; k < vpf; k++)
                {
                    int l, s;
                    uf.find(vpf * 4 * i + vpf * j + k, l, s);
                    if (integerreduceddofs.count(labelmap[l]))
                        intcoords++;
                }
                if (intcoords > maxintcoords)
                {
                    maxintcoords = intcoords;
                    pindof = 4 * vpf * i + vpf * j;
                }
            }
        }
        for (int i = 0; i < vpf; i++)
        {
            int l, s;
            uf.find(pindof + i, l, s);
            Ccoeffs.push_back({ nconstraints, labelmap[l], 1.0 });
            nconstraints++;
        }
        if (verbose)
            std::cout << "Pinning DOFs " << pindof << " through " << pindof + vpf << " to the origin (" << maxintcoords << " integer dimensions)" << std::endl;

        if(verbose)
            std::cout << "Enforcing " << nconstraints << " constraints" << std::endl;
        
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
                    rhs[3 * vpf * i + 3 * l + k] = opt.scale * field.tetFrame(i)(l, k);
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
        
//        if (opt.parameterizationType == CubeCoverOptions::ParameterizationType::PT_SEAMLESS)
        {            
            for (int i = 0; i < jumpdofs; i++)
            {
                intdofs.push_back(reduceddofs + i);
            }
        }
        for (auto it : integerreduceddofs)
        {
            intdofs.push_back(it);
        }

        if (!MIP(C, MIPA, result, MIPrhs, intdofs, opt))
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

            for (int i : intdofs)
                result[i] = std::round(result[i]);

            Eigen::VectorXd augresult(result.size() + 1);
            augresult.segment(0, result.size()) = result;
            augresult[result.size()] = 1;
            std::cout << "C: " << (C * augresult).norm() << std::endl;           

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