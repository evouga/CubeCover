#include "FrameFieldVis.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"

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