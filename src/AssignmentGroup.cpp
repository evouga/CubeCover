#include "AssignmentGroup.h"
#include <utility>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense>

namespace CubeCover {


    AssignmentGroup::AssignmentGroup(int vectorsPerFrame) : vpf(vectorsPerFrame)
    {
        permutation.resize(vectorsPerFrame);
        sign.resize(vectorsPerFrame);
        for (int i = 0; i < vectorsPerFrame; i++)
        {
            permutation[i] = i;
            sign[i] = 1;
        }
    }

    bool AssignmentGroup::isValid() const
    {
        for (int i = 0; i < vectorsPerFrame(); i++)
        {
            bool found = false;
            for (int j = 0; j < vectorsPerFrame(); j++)
            {
                if (permutation[j] == i)
                    found = true;
            }
            if (!found)
                return false;

            if (sign[i] != 1 && sign[i] != -1)
                return false;
        }

        return true;
    }

    bool AssignmentGroup::isIdentity() const
    {
        for (int i = 0; i < vectorsPerFrame(); i++)
        {
            if (permutation[i] != i)
                return false;
            if (sign[i] != 1)
                return false;
        }
        return true;
    }

    AssignmentGroup AssignmentGroup::inverse() const
    {
        AssignmentGroup result(vectorsPerFrame());
        for (int i = 0; i < vectorsPerFrame(); i++)
        {
            result.permutation[permutation[i]] = i;
            result.sign[permutation[i]] = sign[i];
        }
        return result;
    }

    void AssignmentGroup::generateAllElements(int vectorsPerFrame, std::vector<AssignmentGroup>& allelements)
    {
        allelements.clear();
        std::vector<int> perm(vectorsPerFrame);
        for (int i = 0; i < vectorsPerFrame; i++)
            perm[i] = i;

        do {
            AssignmentGroup o(vectorsPerFrame);
            for (int i = 0; i < vectorsPerFrame; i++)
            {
                o.permutation[i] = perm[i];
            }
            for (uint64_t j = 0; j < (1LL << vectorsPerFrame); j++)
            {
                AssignmentGroup o2 = o;
                for (int k = 0; k < vectorsPerFrame; k++)
                    o2.sign[k] = (j & (1LL << k)) ? -1 : 1;
                allelements.push_back(o2);
            }

        } while (std::next_permutation(perm.begin(), perm.end()));
    }

    int AssignmentGroup::orientation() const
    {
        Eigen::MatrixXd P(vectorsPerFrame(),vectorsPerFrame());
        P.setZero();
        for (int i = 0; i < vectorsPerFrame(); i++)
        {
            P(targetVector(i), i) = targetSign(i);
        }
        return P.determinant();
    }

    AssignmentGroup operator*(const AssignmentGroup& second, const AssignmentGroup& first)
    {
        int vpf1 = first.vectorsPerFrame();
        int vpf2 = second.vectorsPerFrame();
        assert(vpf1 == vpf2);
        AssignmentGroup o(vpf1);
        for (int i = 0; i < vpf1; i++)
        {
            o.permutation[i] = second.permutation[first.permutation[i]];
            o.sign[i] = second.sign[first.permutation[i]] * first.sign[i];
        }
        return o;
    }

    Eigen::MatrixXd operator*(const AssignmentGroup& o, const Eigen::MatrixXd& frame)
    {
        int vpf = o.vectorsPerFrame();
        assert(frame.rows() == vpf);
        assert(frame.cols() == 3);
        Eigen::MatrixXd result(vpf, 3);
        for (int i = 0; i < vpf; i++)
        {
            result.row(o.permutation[i]) = frame.row(i) * o.sign[i];
        }
        return result;
    }

    AssignmentGroup localAssignment(const Eigen::MatrixXd& from, const Eigen::MatrixXd& to)
    {
        int vpf = from.rows();
        assert(to.rows() == vpf);
        assert(from.cols() == 3);
        assert(to.cols() == 3);
        AssignmentGroup besto(vpf);
        double bestdist = std::numeric_limits<double>::infinity();
        std::vector<AssignmentGroup> elements;
        AssignmentGroup::generateAllElements(vpf, elements);
        for (auto& it : elements)
        {
            Eigen::MatrixXd mapped = it * from;
            double dist = 0;
            for (int i = 0; i < vpf; i++)
            {
                dist += (mapped.row(i) - to.row(i)).squaredNorm();
            }
            if (dist < bestdist)
            {
                bestdist = dist;
                besto = it;
            }
        }
        return besto;
    }

    std::ostream& operator<<(std::ostream& os, const AssignmentGroup& o)
    {
        os << '[';
        for (int i = 0; i < o.vectorsPerFrame(); i++)
        {
            if (i != 0)
                os << ' ';
            os << o.sign[i] * (1 + o.permutation[i]);
        }
        os << ']';
        return os;
    }

    bool parseFromPerm(const std::vector<int>& pvals, AssignmentGroup& result)
    {
        int vpf = pvals.size();
        AssignmentGroup g(vpf);
        for (int i = 0; i < vpf; i++)
        {
            g.permutation[i] = std::abs(pvals[i]) - 1;
            g.sign[i] = pvals[i] > 0 ? 1 : -1;
        }
        if (!g.isValid())
            return false;
        result = g;
        return true;
    }
};