#include "SingularCurveNetwork.h"
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include <map>
#include "AssignmentGroup.h"
#include <iostream>

void extractSingularCurveNetwork(const Eigen::MatrixXd& V, 
    const CubeCover::TetMeshConnectivity& mesh, 
    const CubeCover::FrameField& field, 
    Eigen::MatrixXd& Pgreen, Eigen::MatrixXi& Egreen,
    Eigen::MatrixXd& Pblue, Eigen::MatrixXi& Eblue,
    Eigen::MatrixXd& Pblack, Eigen::MatrixXi& Eblack
)
{
    std::map<int, int> blacksingverts;
    std::map<int, int> greensingverts;
    std::map<int, int> bluesingverts;
    int nextblackidx = 0;
    int nextblueidx = 0;
    int nextgreenidx = 0;
    std::vector<int> blacksingedges;
    std::vector<int> bluesingedges;
    std::vector<int> greensingedges;
    int nsingedges = field.nSingularEdges();
    int vpf = field.vectorsPerFrame();

    for (int i = 0; i < nsingedges; i++)
    {
        int edge = field.singularEdge(i);
        int v0 = mesh.edgeVertex(edge, 0);
        int v1 = mesh.edgeVertex(edge, 1);

        bool isgreen = false;
        bool isblue = false;

        if (vpf == 3)
        {
            int nbtet = mesh.nEdgeTets(edge);
            CubeCover::AssignmentGroup o(vpf);
            for (int j = 0; j < nbtet; j++)
            {
                int tet = mesh.edgeTet(edge, j);
                int faceidx = mesh.edgeTetFaceIndex(edge, j, 1);
                int face = mesh.tetFace(tet, faceidx);
                int orient = mesh.tetFaceOrientation(tet, faceidx);
                CubeCover::AssignmentGroup faceo = field.faceAssignment(face);

                if (orient == 1)
                    faceo = faceo.inverse();
                o = faceo * o;
            }

            int fixed = -1;
            for (int j = 0; j < 3; j++)
            {
                if (o.targetVector(j) == j && o.targetSign(j) == 1)
                {
                    fixed = j;
                }
            }
            if (fixed != -1 && nbtet > 0)
            {
                Eigen::Vector3d edgevec = V.row(mesh.edgeVertex(edge, 1)).transpose() - V.row(mesh.edgeVertex(edge, 0)).transpose();

                CubeCover::AssignmentGroup o(vpf);
                int votesplus = 0;
                int votesminus = 0;
                for (int j = 0; j < nbtet; j++)
                {
                    int tet = mesh.edgeTet(edge, j);
                    int faceidx = mesh.edgeTetFaceIndex(edge, j, 1);
                    int face = mesh.tetFace(tet, faceidx);
                    int orient = mesh.tetFaceOrientation(tet, faceidx);

                    Eigen::Vector3d framevec = field.tetFrame(tet).row(o.targetVector(fixed)).transpose();
                    framevec *= o.targetSign(fixed);

                    if (framevec.dot(edgevec) > 0)
                    {
                        votesplus++;
                    }
                    else
                    {
                        votesminus++;
                    }

                    CubeCover::AssignmentGroup faceo = field.faceAssignment(face);

                    if (orient == 1)
                        faceo = faceo.inverse();
                    o = faceo * o;
                }
                int sign = (votesplus >= votesminus ? 1 : -1);
                int f1 = (fixed + 1) % 3;
                int f2 = (fixed + 2) % 3;
                if (o.targetVector(f1) == f2 && o.targetSign(f1) == -1 && o.targetVector(f2) == f1 && o.targetSign(f2) == 1)
                {
                    if (sign == 1)
                        isgreen = true;
                    else
                        isblue = true;
                }
                else if (o.targetVector(f1) == f2 && o.targetSign(f1) == 1 && o.targetVector(f2) == f1 && o.targetSign(f2) == -1)
                {
                    if (sign == 1)
                        isblue = true;
                    else
                        isgreen = true;
                }
            }
        }
        if (isgreen)
        {
            auto it = greensingverts.find(v0);
            if (it == greensingverts.end())
            {
                greensingverts[v0] = nextgreenidx++;
            }
            it = greensingverts.find(v1);
            if (it == greensingverts.end())
            {
                greensingverts[v1] = nextgreenidx++;
            }
            greensingedges.push_back(edge);
        }
        else if (isblue)
        {
            auto it = bluesingverts.find(v0);
            if (it == bluesingverts.end())
            {
                bluesingverts[v0] = nextblueidx++;
            }
            it = bluesingverts.find(v1);
            if (it == bluesingverts.end())
            {
                bluesingverts[v1] = nextblueidx++;
            }
            bluesingedges.push_back(edge);
        }
        else
        {
            auto it = blacksingverts.find(v0);
            if (it == blacksingverts.end())
            {
                blacksingverts[v0] = nextblackidx++;
            }
            it = blacksingverts.find(v1);
            if (it == blacksingverts.end())
            {
                blacksingverts[v1] = nextblackidx++;
            }
            blacksingedges.push_back(edge);
        }
    }

    Pblack.resize(blacksingverts.size(), 3);
    for (auto it : blacksingverts)
    {
        Pblack.row(it.second) = V.row(it.first);
    }

    int nblacksingedges = blacksingedges.size();
    Eblack.resize(nblacksingedges, 2);
    for (int i = 0; i < nblacksingedges; i++)
    {
        int edge = blacksingedges[i];
        int v0 = mesh.edgeVertex(edge, 0);
        int v1 = mesh.edgeVertex(edge, 1);
        Eblack(i, 0) = blacksingverts[v0];
        Eblack(i, 1) = blacksingverts[v1];
    }


    Pgreen.resize(greensingverts.size(), 3);
    for (auto it : greensingverts)
    {
        Pgreen.row(it.second) = V.row(it.first);
    }

    int ngreensingedges = greensingedges.size();
    Egreen.resize(ngreensingedges, 2);
    for (int i = 0; i < ngreensingedges; i++)
    {
        int edge = greensingedges[i];
        int v0 = mesh.edgeVertex(edge, 0);
        int v1 = mesh.edgeVertex(edge, 1);
        Egreen(i, 0) = greensingverts[v0];
        Egreen(i, 1) = greensingverts[v1];
    }


    Pblue.resize(bluesingverts.size(), 3);
    for (auto it : bluesingverts)
    {
        Pblue.row(it.second) = V.row(it.first);
    }

    int nbluesingedges = bluesingedges.size();
    Eblue.resize(nbluesingedges, 2);
    for (int i = 0; i < nbluesingedges; i++)
    {
        int edge = bluesingedges[i];
        int v0 = mesh.edgeVertex(edge, 0);
        int v1 = mesh.edgeVertex(edge, 1);
        Eblue(i, 0) = bluesingverts[v0];
        Eblue(i, 1) = bluesingverts[v1];
    }
}