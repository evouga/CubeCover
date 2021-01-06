#ifndef CUBECOVER_H
#define CUBECOVER_H

#include <Eigen/Core>
#include <vector>

namespace CubeCover
{

    struct CubeCoverOptions
    {
        CubeCoverOptions() :
            parameterizationType(ParameterizationType::PT_SEAMLESS),
            assignmentHandling(AssignmentHandling::AH_USEPROVIDED),
            boundaryConditions(BoundaryConditions::BC_FREE),
            MIPtol(1e-6),
            scale(1.0),
            verbose(false)
        {}

        /*
         * Whether to compute a seamless parameterization or integer-grid
         * parameterization. In both cases, the integer-value isosurfaces of the
         * parameterization is continuous across tetrahedral faces. An integer-
         * grid parameterization in addition will have all singular vertices and
         * edges aligned to the integer parameter grid. (This additional property
         * is necessary for meshability near singularities and is required by
         * downstream tools such as HexEx).
         * PT_DISCONTINUOUS simply integrates the frame field on the input
         * tetrahedral mesh, cut into a simply-connected region using an
         * arbitrary set of seam faces (without enforcing continuity of the
         * integer parameter grid across the seam). Useful mainly for 
         * debugging.
         */
        enum struct ParameterizationType
        {
            PT_SEAMLESS,
            PT_INTEGERGRID,
            PT_DISCONTINUOUS
        } parameterizationType;

        /*
         * How singularities in the input frame field are specified.
         * AH_USEPROVIDED instructs CubeCover that local assignments between
         * neighboring frames (across tetrahedral mesh interior faces) have already
         * been computed in advance and that it should respect those assignments.
         * AH_RECOMPUTE asks CubeCover to ignore any provided assignments and to
         * compute assignments itself. It will do so by fitting the best possible
         * frame permutation and sign assignment across each interior face (making
         * no attempt to guarantee that the resulting set of local assignments is
         * topologically valid).
         * Use AH_RECOMPUTE if you only have frames and no priors about the
         * singular structure. But do not expect good results if there is
         * significant noise in the input frame field.
         */
        enum struct AssignmentHandling
        {
            AH_USEPROVIDED,
            AH_RECOMPUTE
        } assignmentHandling;

        /*
         * The boundary conditions to use during parameterization.
         * Free boundaries are unconstrained in any way.
         * Use BC_FORCEINTEGER to require the parameterization to align all
         * boundary faces of the input tetrahedral mesh to the integer grid of
         * the parameter domain. (Obviously, doing so will increase the distortion
         * of the computed parameterization). CubeCover will assign each boundary
         * face of the input mesh to one of the parameter axis based on the closest
         * match between the face normal and the frame field vectors on the face's
         * tetrahedron. It will make no attempt to guarantee a topologically valid
         * or geometrically smooth normal assignment if the input frame field has
         * high noise.
         */
        enum struct BoundaryConditions
        {
            BC_FREE,
            BC_FORCEINTEGER
        } boundaryConditions;

        /*
         * Convergence tolerance passed to the MIP solver.
         */
        double MIPtol;

        /*
         * A unform constant that all frame fields will be scaled by before
         * parameterization. Adjust this parameter to tune the scale (number of
         * integer isosurfaces) in the parameterization.
         */
        double scale;

        bool verbose;
    };

    /*
     * Computes a seamless or integer-grid parameterization given a tetrahedral 
     * mesh and a frame field (discretized on the cells of the tet mesh).
     * Inputs:
     *  V           The nVerts x 3 matrix of vertex positions. 
     *  T           The nTets x 4 matrix of tetrahedra. Each row of T specifies 
     *              the four vertex indices that make up the tet. Vertices are 
     *              zero-indexed. The tetrahedral mesh must be face-connected
     *              and manifold.
     *  frames      k*nTets x 3 matrix specifying the frame field to integrate.
     *              k is the number of vectors that make up each frame (you
     *              will typically want this to be 3).
     *              The first k lines describe the frame field on the first 
     *              tetrahedron in the .mesh file, the second group of k lines 
     *              describes the frame on the second tet, etc.
     *              Each frame is encoded in the order:
     *                  v1x v1y v1z
     *                  v2x v2y v2z
     *                  v3x v3y v3z
     *                  ...
     *              where v1, v2, v3, ... are the vectors (in Euclidean ambient 
     *              coordinates) of the frame.
     *  assignments Information about how the frames on adjacent tetrahedra
     *              relate to each other. Each row encodes a frame assignment 
     *              across one face of the tet mesh. assignments must have
     *              2+k columns, with rows:
     *                  t1 t2 p1 p2 p3 ...
     *              where t1 and t2 are the (zero-indexed) indices of two 
     *              tetrahedra. They must satisfy 0 <= t1,t2 < n and the two 
     *              tetrahedra must be face neighbors in T.
     *              The set {t1, t2} may appear at most once (in either order)
     *              in the assignments matrix. (The assignment coupling two 
     *              neighboring tetrahedra is assumed to be the identity unless 
     *              otherwise specified).
     *              The integers p1, p2, p3, ... must be in the set 
     *                  {-1,-2,-3,...,-k,1,2,3,...,k} 
     *              and they encode a frame assignment as follows: 
     *              vector i (for i = {1,2,3,...,k}) in the frame on 
     *              tetrahedron t1, when multiplied by sign(pi), matches vector 
     *              |pi| on tetrahedron t2. 
     *  opt         Various options that control the parameterization. See above.
     * Output:
     *  parameterization:   An 4*nTets x k matrix of parameter values. Row
     *                      4*i+j lists the parameter values for vertex T(i,j)
     *                      on tetrahedron i. Note that the parameterization is
     *                      defined on a "tetrahedral soup" of the original
     *                      mesh, i.e. the parameter value at a vertex will
     *                      differ across each tetrahedron that shares that 
     *                      vertex. Nevertheless, the parameterization is
     *                      seamless, in the sense that the set of integer
     *                      isosurfaces of each parameter function meet
     *                      along coinciding lines at the faces between tets.
     * Returns whether the parameterization succeeded. More information is
     * might be available on failure if verbose is set to true in opt.
     */
    bool cubeCover(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& T,
        const Eigen::MatrixXd& frames,
        const Eigen::MatrixXi& assignments,
        Eigen::MatrixXd &parameterization,
        CubeCoverOptions opt
    );

    /*
     * As above, but use default parameters (integer-grid parameterization, force
     * integer boundary alignment, and use provided frame assignments).
     */
    bool cubeCover(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& T,
        const Eigen::MatrixXd& frames,
        const Eigen::MatrixXi& assignments,
        Eigen::MatrixXd& parameterization
        );

};

#endif