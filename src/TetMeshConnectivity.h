#ifndef TETMESHCONNECTIVITY_H
#define TETMESHCONNECTIVITY_H

#include <Eigen/Core>
#include <vector>

/*
 * Data structure for computing and manipulating tet mesh connectivity. Takes
 * as input a "libigl-style" tetrahedral cell matrix T, of size ntets x 4.
 * Each row of T lists the vertex indices of one tetrahedral cell.
 */

namespace CubeCover {

    class TetMeshConnectivity
    {
    public:
        TetMeshConnectivity(const Eigen::MatrixXi& T);

        int nTets() const { return T.rows(); }
        int nFaces() const { return F.rows(); }
        int nEdges() const { return E.rows(); }

        /*
         * Maps a local index on a tet (0 <= idx <= 3) to the vertex's global
         * index.
         */
        int tetVertex(int tet, int idx) const { return T(tet, idx); }

        /*
        * Maps a local index on a face (0 <= idx <= 2) to the vertex's global
        * index.
        */
        int faceVertex(int face, int idx) const { return F(face, idx); }

        /*
        * Maps a local index on an edge (0 <= idx <= 1) to the vertex's global
        * index.
        */
        int edgeVertex(int edge, int idx) const { return E(edge, idx); }

        /*
        * Given a face and the local index of a neighboring tetrahedron (tidx
        * \in {0,1}), as well as face local vertex index (0 <= vidx <= 2),
        * returns the local index (in [0,3]) of that vertex on that tet.
        */
        int faceTetVertexIndex(int face, int tidx, int vidx) const { return faceTetVertIndices(2 * face + tidx, vidx); }

        /*
         * Maps a local index on a face (0 <= idx <= 1) to the tet's global
         * index. Returns -1 if there is no such tet (which can only happen 
         * when the face is on the boundary).
         */
        int faceTet(int face, int idx) const { return faceTets(face, idx); }

        /*
         * Maps a local index on a tet (0 <= idx <= 3) to the face's global
         * index.
         */
        int tetFace(int tet, int idx) const { return tetFaces(tet, idx); }

        /*
         * Given a local face index on a tet (0 <= idx <= 3), computes the
         * orientation of that face with respect to the tet.
         * Satisfies the invariant: 
         *  faceTet(tetFace(i,j), tetFaceOrientation(i,j)) == i
         */
        int tetFaceOrientation(int tet, int idx) const { return tetFaceOrientations(tet, idx); }

        /*
         * Given the local index of a vertex on a tet (0 <= idx <= 3), computes
         * the global index of the tetrahedron that borders tet on the face 
         * opposite the given vertex. Returns -1 if there is no such tet
         * (because the opposite face is a boundary face).
         */
        int tetOppositeVertex(int tet, int idx) const;

        /*
         * Maps a local index (0 <= idx <= 5) to the edge's global index.
         */
        int tetEdge(int tet, int idx) const { return tetEdges(tet, idx); }

        /*
         * Returns the number of tetrahedra that share a given edge.
         */
        int nEdgeTets(int edge) const { return edgeTets[edge].size(); }

        /*
         * Returns the global tet index of one tet that neighbors the given 
         * edge. (0 <= idx < nEdgeTets(edge)). The tetrahedra are ordered so
         * that they circulate in the direction of the edge's orientation; if 
         * the edge is a boundary edge, the "gap" occurs before tet 0 (and
         * after tet nEdgeTets(edge)-1).
         */
        int edgeTet(int edge, int idx) const { return edgeTets[edge][idx]; }

        /*
         * Given the local index of a tet neighboring an edge 
         * (0 <= idx <= nEdgeTets(edge)), compute the local index on that tet 
         * of the edge. (The value returned will lie in [0, 5]).
         * Satisfies the invariant:
         *  tetEdge(edgeTet(i,j), indexOnEdgeTets(i,j)) == i
         */
        int edgeIndexOnTet(int edge, int idx) const { return indexOnEdgeTets[edge][idx]; }

        /*
         * Given a local index of an edge's neighboring tet (with
         * 0 <= idx <= nEdgeTets(edge)), two faces of that tet will share the edge.
         * This function returns the local index on the tet of those two faces.
         * fidx specifies which of the two faces: 0 is the "from" face and 1 is the
         * "to" face (so that fidx=1 is the face between tets edgeTet(edge, tidx) and
         * edgeTet(edge, tidx+1), assuming both these tets exist.)
         */
        int edgeTetFaceIndex(int edge, int tidx, int fidx) const { return edgeTetFaceIndices[edge][tidx][fidx]; }

        /*
         * Maps a local index on a face to the global edge index. 
         * (0 <= idx <= 2). Edge i is opposite vertex i on each triangle.
         */
        int faceEdge(int face, int idx) const { return faceEdges(face, idx); }

        bool isBoundaryFace(int face) const;
        bool isBoundaryEdge(int edge) const;

        /* 
         * Checks whether the tet mesh is face- and edge-manifold. Set verbose 
         * to true for more detailed information about any problems that are
         * detected with the tet mesh.
         */
        bool isManifold(bool verbose) const;

        /*
         * Checks whether there is a path from any tetrahedral cell to any 
         * other, passing through cells and the faces between them.
         */
        bool isFaceConnected() const;

    private:
        // Primary data structures
        Eigen::MatrixXi T;
        Eigen::MatrixXi F;
        Eigen::MatrixXi E;

        // Tet-Face connectivity
        // faceTets(tetFaces(i,j), tetFaceOrientations(i,j)) == i
        Eigen::MatrixXi faceTets;
        Eigen::MatrixXi tetFaces;
        Eigen::MatrixXi tetFaceOrientations;
        Eigen::MatrixXi faceTetVertIndices;

        // Tet-Edge connectivity
        // tetEdges(edgeTets[i][j], indexOnEdgeTets[i][j]) = i
        Eigen::MatrixXi tetEdges;
        std::vector<std::vector<int> > edgeTets;
        std::vector<std::vector<Eigen::Vector2i > > edgeTetFaceIndices; // 0 has the same orientation, 1 has the opposite orientation
        std::vector<std::vector<int> > indexOnEdgeTets;

        // Face-Edge connectivity
        // faceEdgeOrientations = 1 if the edge orientation matches the face orientation. -1 otherwise.
        Eigen::MatrixXi faceEdges;
        Eigen::MatrixXi faceEdgeOrientations;

    };

};

#endif