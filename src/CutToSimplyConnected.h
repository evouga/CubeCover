#ifndef CUTTOSPHERE_H
#define CUTTOSPHERE_H

#include <vector>

namespace CubeCover
{

    class TetMeshConnectivity;

    /*
     * Computes the set of faces (seams) along which a tet mesh needs to be cut to
     * make the mesh simply-connected. During this process, the protectedEdges
     * will be treated as if they've been removed from the mesh volume (forming
     * one-dimensional infinitesimal tunnels through the mesh). The ouput faces
     * will form a 2-complex whose boundary edges either lie on the boundary of the
     * original tet mesh, or are in protectedEdges.
     *
     * Note that cutting along the cutFaces does not always yield a topological
     * sphere (since this function does not try to "bore out" interior voids).
     *
     * Nothing is guaranteed regarding optimality of the cutFaces.
     *
     * The input tet mesh must be manifold and face-connected.
     */

    enum struct CST_RESULT
    {
        CST_OK = 0,
        CST_NOTMANIFOLD,        // Input mesh is not manfiold. You might get more
                                // details about the problem with verbose=true.
        CST_NOTCONNECTED,       // Input mesh is not face-connected.
        CST_BADEDGEIDX          // protectedEdges contain an invalid edge index.
    };

    CST_RESULT cutToSimplyConnected(
        const TetMeshConnectivity& mesh,
        const std::vector<int>& protectedEdges,
        std::vector<int>& cutFaces,
        bool verbose = false);

};

#endif