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
        // check that fraFilename ends in .fra 
        // if ends .fra call v1
        // if ends in .bfra call v2 
        // else return false

        std::string fraExt = fraFilename.substr(fraFilename.find_last_of(".") + 1);
        if (fraExt == "fra")
        {
            return readFrameField_v1(fraFilename, permFilename, T, frames, assignments, verbose);
        }
        else if (fraExt == "bfra")
        {
            return readFrameField_v2(fraFilename, permFilename, T, frames, assignments, verbose);
        }
        else
        {
            if (verbose)
                std::cerr << "Unknown file extension: " << fraExt << std::endl;
            return false;
        }

    }


    bool readFrameField_v1(const std::string& fraFilename, const std::string& permFilename, const Eigen::MatrixXi& T,
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



    bool readFrameField_v2(const std::string& fraFilename, const std::string& permFilename, const Eigen::MatrixXi& T,
        Eigen::MatrixXd& frames,
        Eigen::MatrixXi& assignments,
        bool verbose)
    {
        TetMeshConnectivity tetMesh(T);


        try {
            std::ifstream inFile(fraFilename, std::ios::binary);
            if (!inFile.is_open()) {
                std::cerr << "Error: Unable to open file for reading: " << fraFilename << std::endl;
                return false;
            }

            char* filename = new char[5];
            inFile.read(reinterpret_cast<char*>(&filename), sizeof("FRA 2") ); 

          // check that it does start with FRA 2
            // if (filename != "FRA 2")
            // {
            //     if (verbose)
            //         std::cerr << "Expected: FRA 2, read: " << filename << std::endl;
            //     return false;
            // }


            // Read matrix rows and cols
            int rows, cols, vpe;
            inFile.read(reinterpret_cast<char*>(&rows), sizeof(int));
            inFile.read(reinterpret_cast<char*>(&cols), sizeof(int));
            inFile.read(reinterpret_cast<char*>(&vpe), sizeof(int));

            // Read matrix data
            frames.resize(rows, cols);
            inFile.read(reinterpret_cast<char*>(frames.data()), rows * cols * sizeof(double));
            inFile.close();
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error: Unable to deserialize matrix: " << e.what() << std::endl;
            return false;
        }


    }























    //     std::ifstream ifs(fraFilename);
    //     if (!ifs)
    //     {
    //         if (verbose)
    //             std::cerr << "Cannot open file: " << fraFilename << std::endl;
    //         return false;
    //     }
    //     std::string header;
    //     ifs >> header;
    //     if (header != "FRA")
    //     {
    //         if(verbose)
    //             std::cerr << "Expected: FRA, read: " << header << std::endl;
    //         return false;
    //     }
    //     int version;
    //     ifs >> version;
    //     if (version != 1)
    //     {
    //         if(verbose)
    //             std::cerr << "Only version 1 of .fra files supported (read: " << version << ")" << std::endl;
    //         return false;
    //     }

    //     int nframes, vpt, type;
    //     ifs >> nframes >> vpt >> type;
    //     if (nframes != tetMesh.nTets())
    //     {
    //         if(verbose)
    //             std::cerr << "Mismatch between .fra file and tet mesh: " << nframes << " frames != " << tetMesh.nTets() << " tetrahedra" << std::endl;
    //         return false;
    //     }

    //     if (vpt <= 0)
    //     {
    //         if (verbose)
    //             std::cerr << "Must have at least one vector per tet" << std::endl;
    //         return false;
    //     }

    //     int ntets = tetMesh.nTets();

    //     frames.resize(ntets * vpt, 3);

    //     bool ok = true;
    //     for (int i = 0; i < ntets * vpt; i++)
    //     {
    //         for (int k = 0; k < 3; k++)
    //         {
    //             ifs >> frames(i, k);
    //             if (!ifs)
    //             {
    //                 if (verbose)
    //                     std::cerr << "Error reading frame data" << std::endl;
    //                 return false;
    //             }
    //         }
    //     }

    //     assignments.resize(0, 2 + vpt);

    //     if (permFilename.length() > 0)
    //     {
    //         std::ifstream permfs(permFilename);
    //         if(!permfs)
    //         {
    //             if (verbose)
    //                 std::cerr << "Cannot open file: " << permFilename << std::endl;
    //             return false;
    //         }
    //         permfs >> header >> version;
    //         if (header != "PERM")
    //         {
    //             if(verbose)
    //                 std::cerr << "Expected: PERM, read: " << header << std::endl;
    //             return false;
    //         }
    //         if (version != 1)
    //         {
    //             if(verbose)
    //                 std::cerr << "Only version 1 of .perm files supported (read: " << version << ")" << std::endl;
    //             return false;
    //         }

    //         int ntet, nf, nv;
    //         permfs >> ntet >> nf >> nv;
    //         if (ntet != tetMesh.nTets())
    //         {
    //             if(verbose)
    //                 std::cerr << "Mismatch between .perm file and tet mesh: " << nframes << " != " << tetMesh.nTets() << " tetrahedra" << std::endl;
    //             return false;
    //         }
    //         if (nf < 0)
    //         {
    //             if (verbose)
    //                 std::cerr << "Bad number of assignment faces specified" << std::endl;
    //             return false;
    //         }
    //         if (nv != vpt)
    //         {
    //             if (verbose)
    //                 std::cerr << "Mismatch between .perm file and .fra file: " << nv << " assignments per face != " << vpt << " vectors per frame" << std::endl;
    //             return false;
    //         }

    //         assignments.resize(nf, 2 + nv);

    //         for (int i = 0; i < nf; i++)
    //         {
    //             for (int j = 0; j < 2 + nv; j++)
    //             {
    //                 permfs >> assignments(i, j);
    //                 if (!permfs)
    //                 {
    //                     if(verbose)
    //                         std::cerr << "Error reading assignment data" << std::endl;
    //                     return false;
    //                 }
    //             }
    //         }
    //     }
    //     return true;
    // }    


};