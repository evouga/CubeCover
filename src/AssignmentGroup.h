#ifndef ASSIGNMENTGROUP_H
#define ASSIGNMENTGROUP_H

#include <vector>
#include <Eigen/Core>
#include <ostream>

namespace CubeCover {

    /*
     * Frame assignment group implementation
     * Generalization of the octahedral group to undirected frames containing
     * arbitrarily many 3-vectors: S_n \times Z_2^n, where n = vectorsPerFrame is
     * the number of vectors in the frame. each element of the group permutes the
     * vectors and possibly also flips some of their signs (without regard to
     * preserving the orientation of the frames).
     */

    class AssignmentGroup
    {
    public:
        /*
         * Generates an identity element.
         */
        AssignmentGroup(int vectorsPerFrame);

        int vectorsPerFrame() const { return vpf; }
        bool isValid() const;
        bool isIdentity() const;
        AssignmentGroup inverse() const;

        /*
         * Whether the assignment is orientation-preserving (+1) or reversing
         * (-1).
         */
        int orientation() const;

        int targetVector(int srcidx) const { return permutation[srcidx]; }
        int targetSign(int srcidx) const { return sign[srcidx]; }

        /*
         * Computes all elements of the group. Danger: very expensive for large n!
         * There are n! 2^n elements!
         */
        static void generateAllElements(int vectorsPerFrame, std::vector<AssignmentGroup>& allelements);

        friend AssignmentGroup operator*(const AssignmentGroup& first, const AssignmentGroup& second);

        /*
         * Computes the group action on frame fields. frame must be a
         * vectorsPerFrame x 3 matrix, with each row of the matrix one vector of
         * the frame, and the matrix returned will also have these dimensions.
         */
        friend Eigen::MatrixXd operator*(const AssignmentGroup& o, const Eigen::MatrixXd& frame);

        friend std::ostream& operator<<(std::ostream& os, const AssignmentGroup& o);

        /*
         * Attempts to create an element from the integer string in one row of
         * the .perm file (see file format specification). Returns false
         * if the pvals are malformed.
         */
        friend bool parseFromPerm(const std::vector<int>& pvals, AssignmentGroup &result);

    private:
        int vpf;
        std::vector<int> permutation;
        std::vector<int> sign;
    };

    /*
     * Finds the element of the assignment group g which minimizes
     * ||g*from - to||^2_F (where the norm is the Frobenius norm). from and to
     * must be matrices of size n x 3 for the same n. Each row of the matrices
     * is one vector in the frame.
     * Danger: very expensive if n is large! n! 2^n possible symmetries will
     * be checked by brute force!
     */
    AssignmentGroup localAssignment(const Eigen::MatrixXd& from, const Eigen::MatrixXd& to);


    AssignmentGroup operator*(const AssignmentGroup& first, const AssignmentGroup& second);
    Eigen::MatrixXd operator*(const AssignmentGroup& o, const Eigen::MatrixXd& frame);
    std::ostream& operator<<(std::ostream& os, const AssignmentGroup& o);
    bool parseFromPerm(const std::vector<int>& pvals, AssignmentGroup &result);
    
};

#endif