#include "SingularCurveNetwork.h"
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include <map>

namespace CubeCover
{

    void extractSingularCurveNetwork(const Eigen::MatrixXd &V, const TetMeshConnectivity &mesh, const FrameField& field, Eigen::MatrixXd& P, Eigen::MatrixXi& E)
    {
        std::map<int, int> singverts;
        int nextidx = 0;
        int nsingedges = field.nSingularEdges();
        for (int i = 0; i < nsingedges; i++)
        {
            int edge = field.singularEdge(i);
            int v0 = mesh.edgeVertex(edge, 0);
            int v1 = mesh.edgeVertex(edge, 1);

            auto it = singverts.find(v0);
            if (it == singverts.end())
            {
                singverts[v0] = nextidx++;
            }
            it = singverts.find(v1);
            if (it == singverts.end())
            {
                singverts[v1] = nextidx++;
            }
        }

        P.resize(singverts.size(), 3);
        for (auto it : singverts)
        {
            P.row(it.second) = V.row(it.first);
        }

        E.resize(nsingedges, 2);
        for (int i = 0; i < nsingedges; i++)
        {
            int edge = field.singularEdge(i);
            int v0 = mesh.edgeVertex(edge, 0);
            int v1 = mesh.edgeVertex(edge, 1);
            E(i, 0) = singverts[v0];
            E(i, 1) = singverts[v1];
        }
    }

};