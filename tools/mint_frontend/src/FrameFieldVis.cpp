#include "FrameFieldVis.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <iostream>
#include <Eigen/Geometry> 

void buildFrameVectors(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    double scale,
    Eigen::MatrixXd& centroids,
    std::vector<Eigen::MatrixXd>& frameVectors   
)
{
    int ntets = mesh.nTets();
    int vpf = field.vectorsPerFrame();
    centroids.resize(ntets, 3);
    frameVectors.resize(2*vpf);
    for (int i = 0; i < 2 * vpf; i++)
    {
        frameVectors[i].resize(ntets, 3);
    }

    for (int i = 0; i < ntets; i++)
    {
        Eigen::Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 4; j++)
        {
            int vert = mesh.tetVertex(i, j);
            centroid += V.row(vert).transpose();            
        }
        centroid /= 4.0;
        centroids.row(i) = centroid.transpose();

        for (int j = 0; j < vpf; j++)
        {
            Eigen::Vector3d vec = field.tetFrame(i).row(j).transpose();
            frameVectors[2 * j].row(i) = scale * vec.transpose();
            frameVectors[2 * j + 1].row(i) = -scale * vec.transpose();
        }
    }
}



void computePerVectorCurl(const Eigen::MatrixXd& V,
    const CubeCover::TetMeshConnectivity& mesh,
    const CubeCover::FrameField& field,
    const std::vector<Eigen::MatrixXd>& frameVectors,   
    Eigen::MatrixXd& splitCurl
)
{
    int ntets = mesh.nTets();
    int nverts = V.rows(); 
    int nfaces = mesh.nFaces();
    int vpf = field.vectorsPerFrame();

    // for (int i = 0; i < ntets; i++)
    // { 


    // }

    splitCurl.resize(ntets, 6);

    // Eigen::MatrixXd& 

    for (int i = 0; i < nfaces; i++)
    { 
        if (!mesh.isBoundaryFace(i) && field.faceAssignment(i).isIdentity())
        {
            int tid0 = mesh.faceTet(i, 0);
            int tid1 = mesh.faceTet(i, 1);

            Eigen::Vector3d v0 = V.row(mesh.faceVertex(i, 0));
            Eigen::Vector3d v1 = V.row(mesh.faceVertex(i, 1));
            Eigen::Vector3d v2 = V.row(mesh.faceVertex(i, 2));
            
            Eigen::Vector3d fn = (v1 - v0).cross(v2 - v0);

            Eigen::Vector3d b0 = fn.cross(v1-v0);
            Eigen::Vector3d b1 = fn.cross(b0);
            b0.normalize();
            b1.normalize();


            std::cout << "sanity check should be 1" << b1.norm() << std::endl;

            for (int j = 0; j < 2*vpf; j++)
            {
                Eigen::Vector3d vf0 = frameVectors[j].row(tid0);
                Eigen::Vector3d vf1 = frameVectors[j].row(tid1);

                double faceCurl = (vf0.dot(b0) - vf1.dot(b0)) + (vf0.dot(b1) - vf1.dot(b1));
                splitCurl(tid0,j) += faceCurl*faceCurl;
                splitCurl(tid1,j) += faceCurl*faceCurl;

            }  
        }

        // for (int j = 0; j < vpf; j++)
        // {

        // }  

    }


    // centroids.resize(ntets, 3);
    // frameVectors.resize(2*vpf);


    // loop over tets 
    // loop over edges 
    // store edge id with sign into a list 




}