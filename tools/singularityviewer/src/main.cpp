#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include <Eigen/Core>
#include "ReadFrameField.h"
#include <iostream>
#include <fstream>
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include "SingularCurveNetwork.h"
#include "readMeshFixed.h"
#include "polyscope/surface_mesh.h"
#include "FrameFieldVis.h"
#include "polyscope/point_cloud.h"
#include <random>
#include "WriteFrameField.h"

int main(int argc, char *argv[])
{
    if (argc != 3 && argc != 4 && argc != 5)
    {
        std::cerr << "Usage: singularityviewer (.mesh file) (.fra file) [bad_verts path] [.perm file]" << std::endl;
        return -1;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    
    std::string meshfile = argv[1];
    std::string frafile = argv[2];
    std::string permfile;

    bool recomputeperms = false;

    if (argc == 5)
        permfile = argv[4];
    else
        recomputeperms = true;

    bool showViz = true;
    std::string badverts;

    if (argc >= 4)
    {
        badverts = argv[3];
        if (badverts != "1")
            showViz = false;
    }

    std::cout << badverts << std::endl;


    Eigen::MatrixXi F;
    if (!CubeCover::readMESH(meshfile, V, T, F))
        return -1;

    Eigen::MatrixXd frames;
    Eigen::MatrixXi assignments;
    if (!CubeCover::readFrameField(frafile, permfile, T, frames, assignments, true))
        return -1;

    CubeCover::TetMeshConnectivity mesh(T);
    
    CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
    if (!field)
        return -1;

    if (recomputeperms)
    {
        std::cout << "No face assignments provided, recomputing: ";
        std::cout.flush();
        field->computeLocalAssignments();
        std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
    }
    field->combAssignments();

    Eigen::MatrixXd Pblack;
    Eigen::MatrixXi Eblack;
    Eigen::MatrixXd Pblue;
    Eigen::MatrixXi Eblue;
    Eigen::MatrixXd Pgreen;
    Eigen::MatrixXi Egreen;
    extractSingularCurveNetwork(V, mesh, *field, Pgreen, Egreen, Pblue, Eblue, Pblack, Eblack);

    Eigen::MatrixXd centroids;
    std::vector<Eigen::MatrixXd> framefieldvecs;
    buildFrameVectors(V, mesh, *field, 1.0, centroids, framefieldvecs);

    // make a mesh out of all of the boundary faces
    int nbdry = 0;
    int nfaces = mesh.nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        if (mesh.isBoundaryFace(i))
            nbdry++;
    }
    Eigen::MatrixXi bdryF(nbdry, 3);
    int curidx = 0;
    for (int i = 0; i < nfaces; i++)
    {
        if (mesh.isBoundaryFace(i))
        {
            for (int j = 0; j < 3; j++)
            {
                bdryF(curidx, j) = mesh.faceVertex(i, j);
                
            }
            // fix triangle orientations
            int tet = mesh.faceTet(i, 0);
            if (tet == -1)
            {
                std::swap(bdryF(curidx, 0), bdryF(curidx, 1));
            }
            curidx++;
        }
    }

    // visualize the seams

    std::vector<int> seamfaces;

    int ninverted = 0;
    int nnontrivial = 0;

    for (int i = 0; i < nfaces; i++)
    {
        if (!field->faceAssignment(i).isIdentity())
        {
            seamfaces.push_back(i);
            nnontrivial++;
        }

        if (field->faceAssignment(i).orientation() == -1)
            ninverted++;
    }

    std::cout << "Non-identity face assignments: " << nnontrivial << std::endl;
    if (ninverted > 0)
    {
        std::cout << "Warning: " << ninverted << " face assignments are orientation-reversing" << std::endl;
    }

    int nseamtris = seamfaces.size();

    Eigen::MatrixXd seamV(3 * nseamtris, 3);
    Eigen::MatrixXi seamF(nseamtris, 3);
    for (int i = 0; i < nseamtris; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            seamF(i, j) = 3 * i + j;
            seamV.row(3 * i + j) = V.row(mesh.faceVertex(seamfaces[i], j));
        }
    }

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    if (showViz)
    {
        polyscope::init();

        auto *tetc = polyscope::registerPointCloud("Centroids", centroids);
        glm::vec3 dotcolor(0.1, 0.1, 0.1);
        tetc->setPointColor(dotcolor);
        tetc->setPointRadius(0.001);
        int vpf = framefieldvecs.size();
        for (int i = 0; i < vpf; i++)
        {
            std::stringstream ss;
            ss << "Frame Vector " << i;
            auto *vf = tetc->addVectorQuantity(ss.str(), framefieldvecs[i]);
            vf->setVectorColor({ dist(rng),dist(rng),dist(rng) });
            vf->setVectorRadius(0.001);
            vf->setEnabled(true);
        }

        auto *green = polyscope::registerCurveNetwork("Singular Curves (+1/4)", Pgreen, Egreen);
        green->setColor({ 0.0,1.0,0.0 });

        auto *blue = polyscope::registerCurveNetwork("Singular Curves (-1/4)", Pblue, Eblue);
        blue->setColor({ 0.0,0.0,1.0 });

        auto *black = polyscope::registerCurveNetwork("Singular Curves (irregular)", Pblack, Eblack);
        black->setColor({ 0.0,0.0,0.0 });

        auto *psMesh = polyscope::registerSurfaceMesh("Boundary Mesh", V, bdryF);
        psMesh->setTransparency(0.2);
        psMesh->setSurfaceColor({ 0.5,0.5,0.0 });

        auto* seammesh = polyscope::registerSurfaceMesh("Seam", seamV, seamF);
        seammesh->setSurfaceColor({ 0.0, 0.0, 0.0 });

        // visualize!
        polyscope::show();
    }
    else{
            std::ofstream ofs(badverts);
            if (!ofs)
            {
                return false;
            }
            int nsing = field->nSingularEdges();
            ofs << "ids " << V.rows() << " " << 2 * nsing << std::endl;
            for (int i = 0; i < nsing; i++)
            {
                int edge = field->singularEdge(i);
                int v0 = mesh.edgeVertex(edge, 0);
                int v1 = mesh.edgeVertex(edge, 1);
                ofs << v0 << std::endl << v1 << std::endl;
            }
        
    }

    delete field;
}