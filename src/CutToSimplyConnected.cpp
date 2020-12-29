#include "CutToSimplyConnected.h"
#include "TetMeshConnectivity.h"
#include <set>
#include <deque>
#include <map>
#include <iostream>

namespace CubeCover
{

    CST_RESULT cutToSimplyConnected(
        const TetMeshConnectivity& mesh,
        const std::vector<int>& protectedEdges,
        std::vector<int>& cutFaces,
        bool verbose)
    {
        cutFaces.clear();

        if (!mesh.isManifold(verbose))
            return CST_RESULT::CST_NOTMANIFOLD;
        if (!mesh.isFaceConnected())
            return CST_RESULT::CST_NOTCONNECTED;

        std::set<int> protedges;
        for (auto it : protectedEdges)
        {
            if (it < 0 || it >= mesh.nEdges())
                return CST_RESULT::CST_BADEDGEIDX;
            protedges.insert(it);
        }

        if (mesh.nTets() == 0)
            return CST_RESULT::CST_OK;

        std::vector<bool> deletedTets(mesh.nTets());
        std::vector<bool> deletedFaces(mesh.nFaces());

        // Build a spanning tree on the dual mesh

        // queue of (face, tet) pairs to explore
        std::deque<std::pair<int, int> > q;
        q.push_back({ -1,0 });
        while (!q.empty())
        {
            auto nextpair = q.front();
            q.pop_front();

            if (deletedTets[nextpair.second])
                continue;

            if (nextpair.first != -1)
                deletedFaces[nextpair.first] = true;
            deletedTets[nextpair.second] = true;

            for (int i = 0; i < 4; i++)
            {
                int nbface = mesh.tetFace(nextpair.second, i);
                int nbtet = mesh.faceTet(nbface, 1 - mesh.tetFaceOrientation(nextpair.second, i));
                if (nbtet != -1 && !deletedTets[nbtet])
                    q.push_back({ nbface, nbtet });
            }
        }

        int nfaces = mesh.nFaces();
        // add boundary edges to the protected edges
        for (int i = 0; i < nfaces; i++)
        {
            if (mesh.isBoundaryFace(i))
            {
                for (int j = 0; j < 3; j++)
                {
                    protedges.insert(mesh.faceEdge(i, j));
                }
            }
        }

        // delete boundary faces from spanning tree
        for (int i = 0; i < nfaces; i++)
        {
            if (mesh.isBoundaryFace(i))
                deletedFaces[i] = true;
        }

        // we now have a topological valid cut. Trim spurious triangles from the cut.

        int nedges = mesh.nEdges();
        std::vector<int> edgecounts(nedges);
        int remaining = 0;
        std::map<int, std::vector<int> > edges2face;
        for (int i = 0; i < nfaces; i++)
        {
            if (!deletedFaces[i])
            {
                for (int j = 0; j < 3; j++)
                {
                    edges2face[mesh.faceEdge(i, j)].push_back(i);
                    edgecounts[mesh.faceEdge(i, j)]++;
                }
                remaining++;
            }
        }
        if (verbose)
        {
            std::cout << "Initial topological cut has " << remaining << " faces" << std::endl;
        }

        std::deque<int> todelete;

        for (int i = 0; i < nedges; i++)
        {
            if (!protedges.count(i) && edgecounts[i] == 1)
                todelete.push_back(i);
        }

        while (!todelete.empty())
        {
            int nextedge = todelete.front();
            todelete.pop_front();
            for (auto adjface : edges2face[nextedge])
            {
                if (!deletedFaces[adjface])
                {
                    deletedFaces[adjface] = true;
                    for (int i = 0; i < 3; i++)
                    {
                        int dfedge = mesh.faceEdge(adjface, i);
                        if (!protedges.count(dfedge))
                        {
                            edgecounts[dfedge]--;
                            if (edgecounts[dfedge] == 1)
                                todelete.push_back(dfedge);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < nfaces; i++)
        {
            if (!deletedFaces[i])
                cutFaces.push_back(i);
        }

        return CST_RESULT::CST_OK;
    }

};