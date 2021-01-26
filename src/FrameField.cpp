#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <Eigen/Sparse>
#include <set>
#include <deque>

namespace CubeCover {

    FrameField::FrameField(const TetMeshConnectivity& tetMesh) : mesh(tetMesh)
    {
        int ntets = mesh.nTets();

        vpf = 0;

        frames.resize(ntets);
        for (int i = 0; i < ntets; i++)
        {
            frames[i].resize(0, 3);
        }

        faceperms.resize(mesh.nFaces(), AssignmentGroup(vpf));
    }

    void FrameField::setNoisyEuclidean(double magnitude)
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distOcta(0, 47);
        std::normal_distribution<double> distFrame;

        std::vector<AssignmentGroup> allels;
        AssignmentGroup::generateAllElements(3, allels);

        int ntets = mesh.nTets();
        frames.resize(ntets);
        for (int i = 0; i < ntets; i++)
        {
            frames[i].resize(3, 3);
        }

        vpf = 3;

        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    frames[i](j, k) = ((j == k) ? 1.0 : 0.0) + magnitude * distFrame(rng);
                }
            }
            frames[i] = allels[distOcta(rng)] * frames[i];
        }

        int nfaces = mesh.nFaces();
        for (int i = 0; i < nfaces; i++)
            faceperms[i] = AssignmentGroup(vpf);

        setSingularEdges();

    }

    void FrameField::computeLocalAssignments()
    {
        int nfaces = mesh.nFaces();
        for (int i = 0; i < nfaces; i++)
        {
            if (mesh.isBoundaryFace(i))
            {
                faceperms[i] = AssignmentGroup(vpf);
            }
            else
            {
                int from = mesh.faceTet(i, 0);
                int to = mesh.faceTet(i, 1);
                faceperms[i] = localAssignment(frames[from], frames[to]);
            }
        }

        setSingularEdges();
    }

    void FrameField::setSingularEdges()
    {
        singularEdges.clear();
        int nedges = mesh.nEdges();
        for (int i = 0; i < nedges; i++)
        {
            if (mesh.isBoundaryEdge(i))
                continue;

            int nbtet = mesh.nEdgeTets(i);
            AssignmentGroup o(vpf);
            for (int j = 0; j < nbtet; j++)
            {
                int tet = mesh.edgeTet(i, j);
                int faceidx = mesh.edgeTetFaceIndex(i, j, 1);
                int face = mesh.tetFace(tet, faceidx);
                int orient = mesh.tetFaceOrientation(tet, faceidx);
                AssignmentGroup faceo = faceperms[face];

                if (orient == 1)
                    faceo = faceo.inverse();
                o = faceo * o;

            }

            if (!o.isIdentity())
            {
                singularEdges.push_back(i);
            }
        }
    }

    FrameField* fromFramesAndAssignments(const TetMeshConnectivity& tetMesh, const Eigen::MatrixXd& frames, const Eigen::MatrixXi& assignments, bool verbose)
    {
        int ntets = tetMesh.nTets();
        if (frames.rows() % ntets != 0)
        {
            if (verbose)
            {
                std::cerr << "Number of rows in frames " << frames.rows() << " is not a multiple of the number of tets " << ntets << std::endl;
            }
            return NULL;
        }
        if (frames.rows() <= 0)
        {
            if (verbose)
            {
                std::cerr << "Must specify at least one vector per tet" << std::endl;
            }
            return NULL;
        }
        if (frames.cols() != 3)
        {
            if(verbose)
                std::cerr << "frames must have three columns" << std::endl;
            return NULL;
        }

        int vpt = frames.rows() / ntets;

        if (assignments.cols() != 2 + vpt && assignments.rows() != 0)
        {
            if (verbose)
            {
                std::cerr << "assignments must have 2 + vectorsPerFrame columns" << std::endl;
            }
            return NULL;
        }

        FrameField* field = new FrameField(tetMesh);
        field->vpf = vpt;

        for (int i = 0; i < ntets; i++)
        {
            field->frames[i].resize(vpt, 3);
            for (int j = 0; j < vpt; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    field->frames[i](j, k) = frames(vpt * i + j, k);
                }
            }            
        }

        int nfaces = tetMesh.nFaces();
        for (int i = 0; i < nfaces; i++)
        {
            field->faceperms[i] = AssignmentGroup(vpt);
        }


        std::map<std::pair<int, int>, int> facemap;
        for (int i = 0; i < nfaces; i++)
        {
            int t0 = tetMesh.faceTet(i, 0);
            int t1 = tetMesh.faceTet(i, 1);
            facemap[{t0, t1}] = i;
        }

        std::set<int> alreadydone;

        int nf = assignments.rows();

        for (int i = 0; i < nf; i++)
        {
            int t0 = assignments(i, 0);
            int t1 = assignments(i, 1);
            
            bool ok = false;
            bool invert = false;
            auto it = facemap.find({ t0,t1 });
            if (it == facemap.end())
            {
                it = facemap.find({ t1,t0 });
                invert = true;
                if (it == facemap.end())
                {
                    if (verbose)
                        std::cerr << "No face between tetrahedra: " << t0 << " and " << t1 << std::endl;
                    delete field;
                    return NULL;
                }
            }
            std::vector<int> pvals(vpt);
            for (int j = 0; j < vpt; j++)
            {
                pvals[j] = assignments(i, 2 + j);
            }
            
            AssignmentGroup g(vpt);
            if (!parseFromPerm(pvals, g))
            {
                if (verbose)
                    std::cerr << "Invalid assignment for tetrahedra: " << t0 << " and " << t1 << std::endl;
                delete field;
                return NULL;
            }
            if (invert)
                g = g.inverse();

            if (alreadydone.count(it->second))
            {
                if (verbose)
                    std::cerr << "Duplicate entry for tetrahedra: " << t0 << " and " << t1 << std::endl;
                delete field;
                return NULL;
            }
            alreadydone.insert(it->second);

            field->faceperms[it->second] = g;
        }

        field->setSingularEdges();

        return field;
    }

    void FrameField::combAssignments()
    {
        int ntets = mesh.nTets();
        if (ntets == 0)
            return;

        int vpf = vectorsPerFrame();

        std::vector<bool> visited(ntets);
        struct Visit
        {
            int nextnode;
            AssignmentGroup g;
        };

        AssignmentGroup identity(vpf);
        std::deque<Visit> q;
        q.push_back({ 0,identity });
        visited[0] = true;

        while (!q.empty())
        {
            Visit next = q.front();
            q.pop_front();
            // propagate the permutation

            frames[next.nextnode] = next.g * frames[next.nextnode];
            for (int i = 0; i < 4; i++)
            {
                int face = mesh.tetFace(next.nextnode, i);
                int orient = mesh.tetFaceOrientation(next.nextnode, i);
                int nb = mesh.faceTet(face, 1 - orient);
                if (nb == -1)
                    continue;

                AssignmentGroup nextg(vpf);
                
                if (orient == 0)
                {
                    faceperms[face] = faceperms[face] * next.g.inverse();                    
                    nextg = faceperms[face].inverse();
                }
                else
                {
                    faceperms[face] = next.g * faceperms[face];
                    nextg = faceperms[face];
                }

                if (!visited[nb])
                {
                    q.push_back({ nb, nextg });
                    visited[nb] = true;
                }
            }
        }
    }
};