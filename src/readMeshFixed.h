#ifndef READMESHFIXED_H
#define READMESHFIXED_H

#include <string>
#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include <iostream>

namespace CubeCover {

    template <typename DerivedV, typename DerivedF, typename DerivedT>
    bool readMESH(
        const std::string mesh_file_name,
        Eigen::PlainObjectBase<DerivedV>& V,
        Eigen::PlainObjectBase<DerivedT>& T,
        Eigen::PlainObjectBase<DerivedF>& F)
    {
        V.resize(0, 0);
        T.resize(0, 0);
        F.resize(0, 0);
        auto nextNonComment = [](std::istream& fs, std::string& s) -> bool
        {
            while (true)
            {
                std::getline(fs, s);
                if (!fs)
                    return false;
                for (int i = 0; i < s.length(); i++)
                {
                    if (s[i] == '#')
                        break;
                    else if (!std::isspace(s[i]))
                        return true;
                }
            }
            // unreachable
            return false;
        };

        std::ifstream ifs(mesh_file_name);
        if (!ifs)
        {
            std::cerr << "Error: couldn't read .mesh file: " << mesh_file_name << std::endl;
            return false;
        }

        std::string line;

        // check that first word is MeshVersionFormatted
        if (!nextNonComment(ifs, line))
        {
            std::cerr << "Error: missing MeshVersionFormatted line" << std::endl;
            return false;
        }
        std::stringstream ss(line);
        std::string format;
        int versionid;
        ss >> format >> versionid;
        if (!ss || format != "MeshVersionFormatted")
        {
            std::cerr << "Error: expected \"MeshVersionFormatted (1 or 2)\", read " << line << std::endl;
            return false;
        }
        if (versionid != 1 && versionid != 2)
        {
            std::cerr << "Error: MeshVersionFormatted must be 1 or 2" << std::endl;
            return false;
        }

        // check that third word is Dimension
        if (!nextNonComment(ifs, line))
        {
            std::cerr << "Error: missing Dimension line" << std::endl;
            return false;
        }
        ss = std::stringstream(line);
        std::string dimension;
        ss >> dimension;
        if (!ss || dimension != "Dimension")
        {
            std::cerr << "Error: expected \"Dimension\", read " << line << std::endl;
            return false;
        }
        int dim;
        ss >> dim;
        if (!ss)
        {
            // not spec compliant, but occurs in the wild
            if (!nextNonComment(ifs, line))
            {
                std::cerr << "Error: missing dimension number" << std::endl;
                return false;
            }
            ss = std::stringstream(line);
            ss >> dim;
            if (!ss)
            {
                std::cerr << "Error: missing dimension number" << std::endl;
                return false;
            }
        }
        if (dim != 3)
        {
            std::cerr << "Error: only Dimension 3 supported not " << dim << std::endl;
            return false;
        }

        while (true)
        {
            if (!nextNonComment(ifs, line))
            {
                std::cerr << "Error: file ended before \"End\" line" << std::endl;
                return false;
            }

            ss = std::stringstream(line);
            std::string sectionid;
            ss >> sectionid;
            if (sectionid == "End")
                return true;
            else if (sectionid == "Vertices")
            {
                // Not spec compliant, but since this was in the IGL code, needed in some cases?
                int nverts;
                ss >> nverts;
                if (!ss)
                {
                    // On next line
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: missing number of vertices following \"Vertices\"" << std::endl;
                        return false;
                    }
                    ss = std::stringstream(line);
                    ss >> nverts;
                    if (!ss)
                    {
                        std::cerr << "Error: missing number of vertices following \"Vertices\"" << std::endl;
                        return false;
                    }
                }
                V.resize(nverts, 3);
                for (int i = 0; i < nverts; i++)
                {
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: read " << i << " vertices, was expecting " << nverts << std::endl;
                        return false;
                    }
                    int ref;
                    ss = std::stringstream(line);
                    ss >> V(i, 0) >> V(i, 1) >> V(i, 2) >> ref;
                    if (!ss)
                    {
                        std::cerr << "Error: expected \"x y z ref\", read " << line << std::endl;
                        return false;
                    }
                }
            }
            else if (sectionid == "Triangles")
            {
                // Not spec compliant, but since this was in the IGL code, needed in some cases?
                int ntris;
                ss >> ntris;
                if (!ss)
                {
                    // On next line
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: missing number of triangles following \"Triangles\"" << std::endl;
                        return false;
                    }
                    ss = std::stringstream(line);
                    ss >> ntris;
                    if (!ss)
                    {
                        std::cerr << "Error: missing number of triangles following \"Triangles\"" << std::endl;
                        return false;
                    }
                }
                F.resize(ntris, 3);
                for (int i = 0; i < ntris; i++)
                {
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: read " << i << " triangles, was expecting " << ntris << std::endl;
                        return false;
                    }
                    int ref;
                    typename DerivedF::Scalar v1, v2, v3;
                    ss = std::stringstream(line);
                    ss >> v1 >> v2 >> v3 >> ref;
                    if (!ss)
                    {
                        std::cerr << "Error: expected \"v1 v2 v3 ref\", read " << line << std::endl;
                        return false;
                    }
                    F(i, 0) = v1 - 1;
                    F(i, 1) = v2 - 1;
                    F(i, 2) = v3 - 1;
                }
            }

            else if (sectionid == "Quadrilaterals")
            {
                // Not spec compliant, but since this was in the IGL code, needed in some cases?
                int nquads;
                ss >> nquads;
                if (!ss)
                {
                    // On next line
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: missing number of quadrilaterals following \"Quadrilaterals\"" << std::endl;
                        return false;
                    }
                    ss = std::stringstream(line);
                    ss >> nquads;
                    if (!ss)
                    {
                        std::cerr << "Error: missing number of quadrilaterals following \"Quadrilaterals\"" << std::endl;
                        return false;
                    }
                }
                if (nquads != 0)
                {
                    std::cerr << "Error: Only tetrahedral meshes are supported (file contains Quadrilaterals)" << std::endl;
                    return false;
                }
            }
            else if (sectionid == "Tetrahedra")
            {
                // Not spec compliant, but since this was in the IGL code, needed in some cases?
                int ntets;
                ss >> ntets;
                if (!ss)
                {
                    // On next line
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: missing number of tetrahedra following \"Tetrahedra\"" << std::endl;
                        return false;
                    }
                    ss = std::stringstream(line);
                    ss >> ntets;
                    if (!ss)
                    {
                        std::cerr << "Error: missing number of tetrahedra following \"Tetrahedra\"" << std::endl;
                        return false;
                    }
                }
                T.resize(ntets, 4);
                for (int i = 0; i < ntets; i++)
                {
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: read " << i << " tetrahedra, was expecting " << ntets << std::endl;
                        return false;
                    }
                    int ref;
                    typename DerivedT::Scalar v1, v2, v3, v4;
                    ss = std::stringstream(line);
                    ss >> v1 >> v2 >> v3 >> v4 >> ref;
                    if (!ss)
                    {
                        std::cerr << "Error: expected \"v1 v2 v3 v4 ref\", read " << line << std::endl;
                        return false;
                    }
                    T(i, 0) = v1 - 1;
                    T(i, 1) = v2 - 1;
                    T(i, 2) = v3 - 1;
                    T(i, 3) = v4 - 1;
                }
            }
            else if (sectionid == "Hexahedra")
            {
                // Not spec compliant, but since this was in the IGL code, needed in some cases?
                int nhexes;
                ss >> nhexes;
                if (!ss)
                {
                    // On next line
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: missing number of hexahedra following \"Hexahedra\"" << std::endl;
                        return false;
                    }
                    ss = std::stringstream(line);
                    ss >> nhexes;
                    if (!ss)
                    {
                        std::cerr << "Error: missing number of hexahedra following \"Hexahedra\"" << std::endl;
                        return false;
                    }
                }
                if (nhexes != 0)
                {
                    std::cerr << "Error: Only tetrahedral meshes are supported (file contains Hexahedra)" << std::endl;
                    return false;
                }
            }
            else if (sectionid == "Edges")
            {
                // Not spec compliant, but since this was in the IGL code, needed in some cases?
                int nedges;
                ss >> nedges;
                if (!ss)
                {
                    // On next line
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: missing number of edges following \"Edges\"" << std::endl;
                        return false;
                    }
                    ss = std::stringstream(line);
                    ss >> nedges;
                    if (!ss)
                    {
                        std::cerr << "Error: missing number of edges following \"Edges\"" << std::endl;
                        return false;
                    }
                }
                for (int i = 0; i < nedges; i++)
                {
                    if (!nextNonComment(ifs, line))
                    {
                        std::cerr << "Error: read " << i << " edges, was expecting " << nedges << std::endl;
                        return false;
                    }
                    int e1, e2, ref;
                    ss = std::stringstream(line);
                    ss >> e1 >> e2 >> ref;
                    if (!ss)
                    {
                        std::cerr << "Error: expected \"e1 e2 ref\", read " << line << std::endl;
                        return false;
                    }
                }
            }
            else
            {
                std::cerr << "Error: unknown section ID " << sectionid << std::endl;
                return false;
            }
        }

        // unreachable
        return false;
    }

};

#endif
