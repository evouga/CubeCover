#include "ReadFrameField.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include <fstream>
#include <iostream>

namespace CubeCover
{

    bool readFrameField(const std::string& fraFilename, const std::string& permFilename, const Eigen::MatrixXi& T,
        Eigen::MatrixXd& frames,
        Eigen::MatrixXi& assignments,
        bool verbose)
    {
        TetMeshConnectivity tetMesh(T);

        std::ifstream ifs(fraFilename);
        if (!ifs)
        {
            if (verbose)
                std::cerr << "Cannot open file: " << fraFilename << std::endl;
            return false;
        }
        std::string header;
        ifs >> header;
        if (header != "FRA")
        {
            if(verbose)
                std::cerr << "Expected: FRA, read: " << header << std::endl;
            return false;
        }
        int version;
        ifs >> version;
        if (version != 1)
        {
            if(verbose)
                std::cerr << "Only version 1 of .fra files supported (read: " << version << ")" << std::endl;
            return false;
        }

        int nframes, vpt, type;
        ifs >> nframes >> vpt >> type;
        if (nframes != tetMesh.nTets())
        {
            if(verbose)
                std::cerr << "Mismatch between .fra file and tet mesh: " << nframes << " frames != " << tetMesh.nTets() << " tetrahedra" << std::endl;
            return false;
        }

        if (vpt <= 0)
        {
            if (verbose)
                std::cerr << "Must have at least one vector per tet" << std::endl;
            return false;
        }

        int ntets = tetMesh.nTets();

        frames.resize(ntets * vpt, 3);

        bool ok = true;
        for (int i = 0; i < ntets * vpt; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                ifs >> frames(i, k);
                if (!ifs)
                {
                    if (verbose)
                        std::cerr << "Error reading frame data" << std::endl;
                    return false;
                }
            }
        }

        assignments.resize(0, 2 + vpt);

        if (permFilename.length() > 0)
        {
            std::ifstream permfs(permFilename);
            if(!permfs)
            {
                if (verbose)
                    std::cerr << "Cannot open file: " << permFilename << std::endl;
                return false;
            }
            permfs >> header >> version;
            if (header != "PERM")
            {
                if(verbose)
                    std::cerr << "Expected: PERM, read: " << header << std::endl;
                return false;
            }
            if (version != 1)
            {
                if(verbose)
                    std::cerr << "Only version 1 of .perm files supported (read: " << version << ")" << std::endl;
                return false;
            }

            int ntet, nf, nv;
            permfs >> ntet >> nf >> nv;
            if (ntet != tetMesh.nTets())
            {
                if(verbose)
                    std::cerr << "Mismatch between .perm file and tet mesh: " << nframes << " != " << tetMesh.nTets() << " tetrahedra" << std::endl;
                return false;
            }
            if (nf < 0)
            {
                if (verbose)
                    std::cerr << "Bad number of assignment faces specified" << std::endl;
                return false;
            }
            if (nv != vpt)
            {
                if (verbose)
                    std::cerr << "Mismatch between .perm file and .fra file: " << nv << " assignments per face != " << vpt << " vectors per frame" << std::endl;
                return false;
            }

            assignments.resize(nf, 2 + nv);

            for (int i = 0; i < nf; i++)
            {
                for (int j = 0; j < 2 + nv; j++)
                {
                    permfs >> assignments(i, j);
                    if (!permfs)
                    {
                        if(verbose)
                            std::cerr << "Error reading assignment data" << std::endl;
                        return false;
                    }
                }
            }
        }
        return true;
    }    
};