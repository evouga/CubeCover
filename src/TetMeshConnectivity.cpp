#include "TetMeshConnectivity.h"
#include <set>
#include <algorithm>
#include <utility>
#include <map>
#include <iostream>
#include <cassert>
#include <deque>

namespace CubeCover
{

    TetMeshConnectivity::TetMeshConnectivity(const Eigen::MatrixXi& T) : T(T)
    {
        struct Triple
        {
            Triple(int a, int b, int c) : first(a), second(b), third(c)
            {
                if (first > second)
                    std::swap(first, second);
                if (first > third)
                    std::swap(first, third);
                if (second > third)
                    std::swap(second, third);
            }

            bool operator<(const Triple& other) const
            {
                if (first < other.first) return true;
                else if (first > other.first) return false;
                else if (second < other.second) return true;
                else if (second > other.second) return false;
                else return third < other.third;
            }

            int first, second, third;
        };

        int ntets = T.rows();
        std::set<Triple> faces;
        std::set<std::pair<int, int> > edges;

        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int jp1 = (j + 1) % 4;
                int jp2 = (j + 2) % 4;
                faces.insert({ T(i,j), T(i,jp1), T(i,jp2) });
                for (int k = j + 1; k < 4; k++)
                {
                    int v0 = T(i, j);
                    int v1 = T(i, k);
                    if (v0 > v1)
                        std::swap(v0, v1);
                    edges.insert({ v0,v1 });
                }
            }
        }

        int nfaces = faces.size();
        F.resize(nfaces, 3);
        faceTets.resize(nfaces, 2);
        faceTets.setConstant(-1);
        tetFaces.resize(ntets, 4);
        tetFaceOrientations.resize(ntets, 4);
        faceTetVertIndices.resize(2 * nfaces, 3);
        faceTetVertIndices.setConstant(-1);

        std::map<Triple, int> facemap;
        int idx = 0;
        for (auto& it : faces)
        {
            F(idx, 0) = it.first;
            F(idx, 1) = it.second;
            F(idx, 2) = it.third;
            facemap[it] = idx;
            idx++;
        }
        assert(idx == nfaces);

        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int jp0 = j;
                int jp1 = (j + 1) % 4;
                int jp2 = (j + 2) % 4;
                int jp3 = (j + 3) % 4;
                int v0 = T(i, j);
                int v1 = T(i, jp1);
                int v2 = T(i, jp2);
                int sign = j % 2;
                if (v0 > v1)
                {
                    std::swap(v0, v1);
                    std::swap(jp0, jp1);
                    sign ^= 1;
                }
                if (v0 > v2)
                {
                    std::swap(v0, v2);
                    std::swap(jp0, jp2);
                    sign ^= 1;
                }
                if (v1 > v2)
                {
                    std::swap(v1, v2);
                    std::swap(jp1, jp2);
                    sign ^= 1;
                }
                int faceid = facemap[{v0, v1, v2}];
                tetFaces(i, jp3) = faceid;
                tetFaceOrientations(i, jp3) = sign;
                faceTets(faceid, sign) = i;
                faceTetVertIndices(2 * faceid + sign, 0) = jp0;
                faceTetVertIndices(2 * faceid + sign, 1) = jp1;
                faceTetVertIndices(2 * faceid + sign, 2) = jp2;
            }
        }

        int nedges = edges.size();
        std::map<std::pair<int, int>, int> edgemap;
        E.resize(nedges, 2);
        idx = 0;
        for (auto& it : edges)
        {
            E(idx, 0) = it.first;
            E(idx, 1) = it.second;
            edgemap[{it.first, it.second}] = idx;
            idx++;
        }
        assert(idx == nedges);

        faceEdges.resize(nfaces, 3);
        faceEdgeOrientations.resize(nfaces, 3);
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int jp1 = (j + 1) % 3;
                int jp2 = (j + 2) % 3;
                int v1 = F(i, jp1);
                int v2 = F(i, jp2);
                int sign = 1;
                if (v1 > v2)
                {
                    std::swap(v1, v2);
                    sign = -1;
                }
                faceEdges(i, j) = edgemap[{v1, v2}];
                faceEdgeOrientations(i, j) = sign;
            }
        }

        tetEdges.resize(ntets, 6);
        edgeTets.resize(nedges);
        indexOnEdgeTets.resize(nedges);
        for (int i = 0; i < ntets; i++)
        {
            int idx = 0;
            for (int j = 0; j < 4; j++)
            {
                for (int k = j + 1; k < 4; k++)
                {
                    int v1 = T(i, j);
                    int v2 = T(i, k);
                    if (v1 > v2)
                        std::swap(v1, v2);
                    int eid = edgemap[{v1, v2}];
                    tetEdges(i, idx) = eid;
                    edgeTets[eid].push_back(i);
                    indexOnEdgeTets[eid].push_back(idx);
                    idx++;
                }
            }
        }

        edgeTetFaceIndices.resize(nedges);
        for (int i = 0; i < nedges; i++)
        {
            int nnbs = edgeTets[i].size();
            edgeTetFaceIndices[i].resize(nnbs);
            int v1 = E(i, 0);
            int v2 = E(i, 1);
            for (int j = 0; j < nnbs; j++)
            {
                int tet = edgeTets[i][j];
                for (int k = 0; k < 4; k++)
                {
                    int face = tetFace(tet, k);
                    int orient = tetFaceOrientation(tet, k);
                    for (int l = 0; l < 3; l++)
                    {
                        int lp1 = (l + 1) % 3;
                        if (F(face, l) == v1 && F(face, lp1) == v2)
                        {
                            int forient = orient;
                            edgeTetFaceIndices[i][j][forient] = k;
                        }
                        else if (F(face, l) == v2 && F(face, lp1) == v1)
                        {
                            int forient = orient ^ 1;
                            edgeTetFaceIndices[i][j][forient] = k;
                        }
                    }
                }
            }
        }

        // sort the fans around each edge
        for (int i = 0; i < nedges; i++)
        {
            std::map<int, int> invtetmap;
            for (int j = 0; j < edgeTets[i].size(); j++)
            {
                invtetmap[edgeTets[i][j]] = j;
            }
            std::set<int> destfaces;
            for (int j = 0; j < edgeTetFaceIndices[i].size(); j++)
            {
                destfaces.insert(tetFaces(edgeTets[i][j], edgeTetFaceIndices[i][j][1]));
            }
            std::vector<int> starts;
            for (int j = 0; j < edgeTetFaceIndices[i].size(); j++)
            {
                if (!destfaces.count(tetFaces(edgeTets[i][j], edgeTetFaceIndices[i][j][0])))
                    starts.push_back(j);
            }
            std::vector<int> newtets;

            if (starts.empty())
                starts.push_back(0);

            for (auto it : starts)
            {
                int tet = edgeTets[i][it];
                newtets.push_back(tet);
                int nextface = edgeTetFaceIndices[i][invtetmap[tet]][1];
                int nexttet = faceTets(tetFaces(tet, nextface), 1 - tetFaceOrientations(tet, nextface));
                while (nexttet != -1 && nexttet != edgeTets[i][it])
                {
                    tet = nexttet;
                    newtets.push_back(tet);
                    nextface = edgeTetFaceIndices[i][invtetmap[tet]][1];
                    nexttet = faceTets(tetFaces(tet, nextface), 1 - tetFaceOrientations(tet, nextface));
                }
            }
            assert(newtets.size() == edgeTets[i].size());
            std::vector<int> newedgetets(newtets.size());
            std::vector<Eigen::Vector2i> newedgetetfaces(newtets.size());
            std::vector<int> newindexonedgetets(newtets.size());
            for (int j = 0; j < newtets.size(); j++)
            {
                newedgetets[j] = newtets[j];
                newedgetetfaces[j] = edgeTetFaceIndices[i][invtetmap[newtets[j]]];
                newindexonedgetets[j] = indexOnEdgeTets[i][invtetmap[newtets[j]]];
            }
            edgeTets[i] = newedgetets;
            edgeTetFaceIndices[i] = newedgetetfaces;
            indexOnEdgeTets[i] = newindexonedgetets;
        }
    }

    int TetMeshConnectivity::tetOppositeVertex(int tet, int idx) const
    {
        return faceTet(tetFace(tet, idx), 1 - tetFaceOrientation(tet, idx));
    }

    bool TetMeshConnectivity::isBoundaryFace(int face) const
    {
        return faceTet(face, 0) == -1 || faceTet(face, 1) == -1;
    }

    bool TetMeshConnectivity::isBoundaryEdge(int edge) const
    {
        int tet0 = edgeTet(edge, 0);
        int prevface = edgeTetFaceIndex(edge, 0, 0);
        int gface = tetFace(tet0, prevface);
        int orient = 1 - tetFaceOrientation(tet0, prevface);
        return faceTet(gface, orient) == -1;
    }

    bool TetMeshConnectivity::isManifold(bool verbose) const
    {
        // check mesh is face-manifold
        std::map<int, int> facecount;
        int ntets = nTets();
        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                facecount[tetFace(i, j)]++;
            }
        }
        for (auto it : facecount)
        {
            if (it.second > 2)
            {
                if (verbose)
                    std::cerr << "Mesh has non-manifold face " << faceVertex(it.first, 0) << ", " << faceVertex(it.first, 1) << ", " << faceVertex(it.first, 2) << std::endl;
                return false;
            }
        }

        int nfaces = nFaces();

        // check face connectivity data structures
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                assert(faceTet(i, j) >= -1 && faceTet(i, j) < ntets);
                if (faceTet(i, j) == -1)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        assert(faceTetVertexIndex(i, j, k) == -1);
                    }
                }
                else
                {
                    for (int k = 0; k < 3; k++)
                    {
                        int vidx = faceTetVertexIndex(i, j, k);
                        assert(vidx >= 0 && vidx < 4);
                        assert(faceVertex(i, k) == tetVertex(faceTet(i, j), vidx));
                    }
                }
            }
        }
        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                assert(tetFace(i, j) >= 0 && tetFace(i, j) < nfaces);
                assert(tetFaceOrientation(i, j) == 0 || tetFaceOrientation(i, j) == 1);
                assert(i == faceTet(tetFace(i, j), tetFaceOrientation(i, j)));
            }
        }

        // check edge connectivity data structures
        int nedges = nEdges();
        for (int i = 0; i < ntets; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                assert(tetEdge(i, j) >= 0 && tetEdge(i, j) < nedges);
            }
        }
        for (int i = 0; i < nedges; i++)
        {
            for (int j = 0; j < nEdgeTets(i); j++)
            {
                int nbtet = edgeTet(i, j);
                assert(nbtet >= 0 && nbtet < ntets);
                int nbidx = edgeIndexOnTet(i, j);
                assert(nbidx >= 0 && nbidx < 6);
                assert(tetEdge(nbtet, nbidx) == i);
                for (int k = 0; k < 2; k++)
                {
                    int etf = edgeTetFaceIndex(i, j, k);
                    assert(etf >= 0 && etf < 4);
                    int face = tetFace(nbtet, etf);
                    bool ok = false;
                    for (int l = 0; l < 3; l++)
                    {
                        if (faceEdge(face, l) == i)
                            ok = true;
                    }
                    assert(ok);
                }
            }
            int components = 0;
            for (int j = 0; j < nEdgeTets(i); j++)
            {
                int jp1 = (j + 1) % nEdgeTets(i);
                int nextface = edgeTetFaceIndex(i, j, 1);
                int nexttet = faceTet(tetFace(edgeTet(i, j), nextface), 1 - tetFaceOrientation(edgeTet(i, j), nextface));
                assert(nexttet == -1 || nexttet == edgeTet(i, jp1));
                if (nexttet == -1)
                    components++;
            }
            if (components > 1)
            {
                assert(isBoundaryEdge(i));
                if (verbose)
                {
                    std::cerr << "Mesh has non-manifold edge " << edgeVertex(i, 0) << " -- " << edgeVertex(i, 1) << std::endl;
                }
                return false;
            }
            else
            {
                assert(isBoundaryEdge(i) ? components == 1 : components == 0);
            }
        }

        return true;
    }

    bool TetMeshConnectivity::isFaceConnected() const
    {
        int ntets = nTets();

        // check mesh is face-connected
        if (ntets > 0)
        {
            std::vector<bool> visited(ntets);
            std::deque<int> q;
            q.push_back(0);
            visited[0] = true;
            while (!q.empty())
            {
                int next = q.front();
                q.pop_front();
                for (int i = 0; i < 4; i++)
                {
                    int nb = tetOppositeVertex(next, i);
                    if (nb != -1 && !visited[nb])
                    {
                        q.push_back(nb);
                        visited[nb] = true;
                    }
                }
            }
            for (int i = 0; i < ntets; i++)
            {
                if (!visited[i])
                {
                    return false;
                }
            }
        }
        return true;
    }

};