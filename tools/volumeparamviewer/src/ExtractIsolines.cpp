#include "ExtractIsolines.h"
#include "TetMeshConnectivity.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"


struct Segment
{
    Eigen::Vector3d end1;
    Eigen::Vector3d end2;
};


void extractPoints(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
    std::vector<glm::vec3>& points,
    std::vector<std::array<double, 3>>& colors)
{

    points.clear();
    colors.clear();

 //   std::cout << V.size() << " rows: " << V.rows()  << " cols: " << V.cols() << std::endl;
    Eigen::MatrixXd vert_vals( V.rows(), V.cols() );


    // generate points
    for (size_t i = 0; i < V.rows(); i++) {
      points.push_back( glm::vec3{V(i,0), V(i,1), V(i,2)} );
      colors.push_back( { V(i,0), V(i,1), V(i,2) } );
    }

    double minvals[3];
    double maxvals[3];
    for (int j = 0; j < 3; j++)
    {
        minvals[j] = std::numeric_limits<double>::infinity();
        maxvals[j] = -std::numeric_limits<double>::infinity();
    }

    int ntets = mesh.nTets();
    for (int i = 0; i < ntets; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                minvals[k] = std::min(minvals[k], values(4 * i + j, k));
                maxvals[k] = std::max(maxvals[k], values(4 * i + j, k));
            }
        }

        int tosamplelower[3];
        int tosampleupper[3];

        for (int j = 0; j < 3; j++)
        {
            tosamplelower[j] = std::ceil(minvals[j]);
            tosampleupper[j] = std::floor(maxvals[j]);
        }

        for ( int tid = 0; tid < 4; tid++ )
        {
            int vid = mesh.tetVertex(i, tid);
            vert_vals.row(vid) = values.row(4 * i + tid); 

            // vert_vals.row(vid) = V.row(vid); 
        }

    }

    
    points.clear();
    colors.clear();

    for (size_t i = 0; i < V.rows(); i++) {
      // colors.push_back( { V(i,0), V(i,1), V(i,2) } );
      // colors.push_back( { vert_vals(i,0), vert_vals(i,1), vert_vals(i,2) } );
      float rescale = 6.28f;
      float v0 = ( vert_vals(i,0) - minvals[0] ) / rescale;
      float v1 = ( vert_vals(i,1) - minvals[1] ) / rescale;
      float v2 = ( vert_vals(i,2) - minvals[2] ) / rescale;
      v0 = v0 - floor(v0);
      v1 = v1 - floor(v1);
      v2 = v2 - floor(v2);

      std::cout << v0 << " " << v1 << " " << v2 << " " << std::endl;

      colors.push_back( { v0, .5, v1 } );
      points.push_back( glm::vec3{ vert_vals(i,0), vert_vals(i,1), vert_vals(i,2) } );


    }

// visualize!

}


void extractIsolines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const Eigen::MatrixXd& values,
    Eigen::MatrixXd& P,
    Eigen::MatrixXi& E,
    Eigen::MatrixXd& P2,
    Eigen::MatrixXi& E2,
    std::vector<int> badverts)
{
    constexpr double clamptol = 1e-6;

    assert(values.rows() == mesh.nTets() * 4);
    assert(values.cols() == 3);

    std::vector<Segment> segs;
    std::vector<Segment> segs2;
    badverts.clear();

    int ntets = mesh.nTets();
    for (int i = 0; i < ntets; i++)
    {
        double minvals[3];
        double maxvals[3];
        for (int j = 0; j < 3; j++)
        {
            minvals[j] = std::numeric_limits<double>::infinity();
            maxvals[j] = -std::numeric_limits<double>::infinity();
        }
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                minvals[k] = std::min(minvals[k], values(4 * i + j, k));
                maxvals[k] = std::max(maxvals[k], values(4 * i + j, k));
            }
        }

        int tosamplelower[3];
        int tosampleupper[3];

        for (int j = 0; j < 3; j++)
        {
            tosamplelower[j] = std::ceil(minvals[j]);
            tosampleupper[j] = std::floor(maxvals[j]);
        }

		Eigen::Vector3d p0 = values.row(4 * i + 0);
		Eigen::Vector3d p1 = values.row(4 * i + 1);
		Eigen::Vector3d p2 = values.row(4 * i + 2);
		Eigen::Vector3d p3 = values.row(4 * i + 3);

		Eigen::Vector3d d1 = p1 - p0;
		Eigen::Vector3d d2 = p2 - p0;
		Eigen::Vector3d d3 = p3 - p0;

		double vol = d1.cross(d2).dot(d3);
		if (std::fabs(vol) < .000001)
		{
            badverts.push_back(i);
            std::cout << i << " " << std::endl;
			for (int j = 0; j < 4; j++)
			{
                // badverts.push_back(mesh.tetVertex(i,j));
				for (int k = j + 1; k < 4; k++)
				{
					Segment s;
					s.end1 = V.row(mesh.tetVertex(i, j));
					s.end2 = V.row(mesh.tetVertex(i, k));
					segs2.push_back(s);
				}
			}
		}

		for (int j = 0; j < 4; j++)
		{
			Eigen::Vector3d p0 = values.row(4 * i + 0);

			
		}

        // add segments for each direction
        for (int dir = 0; dir < 3; dir++)
        {            
            int sample1 = (dir + 1) % 3;
            int sample2 = (dir + 2) % 3;
            for (int val1 = tosamplelower[sample1]; val1 <= tosampleupper[sample1]; val1++)
            {
                for (int val2 = tosamplelower[sample2]; val2 <= tosampleupper[sample2]; val2++)
                {
                    std::vector<Eigen::Vector4d> pts;
                    int nvertpts = 0;
                    int nedgepts = 0;
                    int nfacepts = 0;

                    // vertex case
                    for (int vert = 0; vert < 4; vert++)
                    {
                        double s1val = values(4 * i + vert, sample1);
                        double s2val = values(4 * i + vert, sample2);
                        if (s1val == double(val1) && s2val == double(val2))
                        {
                            Eigen::Vector4d pt(0, 0, 0, 0);
                            pt[vert] = 1.0;
                            pts.push_back(pt);
                            nvertpts++;
                        }
                    }

                    // edge case
                    for (int vert1 = 0; vert1 < 4; vert1++)
                    {
                        for (int vert2 = vert1 + 1; vert2 < 4; vert2++)
                        {
                            double s11val = values(4 * i + vert1, sample1);
                            double s12val = values(4 * i + vert2, sample1);
                            double s21val = values(4 * i + vert1, sample2);
                            double s22val = values(4 * i + vert2, sample2);
                            if (s11val == double(val1) && s12val == double(val1))
                            {
                                // whole edge intersect the isosurface
                                if (s21val != s22val)
                                {
                                    // equality is handled in the vertex case
                                    double bary = (double(val2) - s21val) / (s22val - s21val);
                                    if (bary > 0.0 && bary < 1.0)
                                    {                                        
                                        Eigen::Vector4d pt(0, 0, 0, 0);
                                        pt[vert1] = 1.0 - bary;
                                        pt[vert2] = bary;
                                        pts.push_back(pt);
                                        nedgepts++;
                                    }
                                }
                                else
                                {
                                    double bary1 = (double(val1) - s11val) / (s12val - s11val);
                                    double bary2 = (double(val2) - s21val) / (s22val - s21val);
                                    if (bary1 > 0.0 && bary1 < 1.0 && bary1 == bary2)
                                    {
                                        Eigen::Vector4d pt(0, 0, 0, 0);
                                        pt[vert1] = 1.0 - bary1;
                                        pt[vert2] = bary1;
                                        pts.push_back(pt);
                                        nedgepts++;
                                    }
                                }
                            }
                        }
                    }


                    // face case
                    for (int vert1 = 0; vert1 < 4; vert1++)
                    {
                        for (int vert2 = vert1 + 1; vert2 < 4; vert2++)
                        {
                            for (int vert3 = vert2 + 1; vert3 < 4; vert3++)
                            {
                                double s11val = values(4 * i + vert1, sample1);
                                double s12val = values(4 * i + vert2, sample1);
                                double s13val = values(4 * i + vert3, sample1);
                                double s21val = values(4 * i + vert1, sample2);
                                double s22val = values(4 * i + vert2, sample2);
                                double s23val = values(4 * i + vert3, sample2);
                                
                                Eigen::Matrix2d M;
                                M << s12val - s11val, s13val - s11val,
                                    s22val - s21val, s23val - s21val;

                                if (M.determinant() == 0)
                                {
									std::cout << "zero volume" << std::endl;
                                    continue;
                                }

                                Eigen::Vector2d rhs(double(val1) - s11val, double(val2) - s21val);
                                Eigen::Vector2d barys = M.inverse() * rhs;
                                if (barys[0] >= 0 && barys[1] >= 0 && barys[0] + barys[1] <= 1.0)
                                {
                                    Eigen::Vector4d pt(0, 0, 0, 0);
                                    pt[vert1] = 1.0 - barys[0] - barys[1];
                                    pt[vert2] = barys[0];
                                    pt[vert3] = barys[1];
                                    pts.push_back(pt);
                                    nfacepts++;
                                }                                
                            }
                        }
                    }

                    // clamp together identical pts
                    std::vector<int> tokeep;
                    int ncrossings = pts.size();
                    for (int j = 0; j < ncrossings; j++)
                    {
                        bool ok = true;
                        for (int k = 0; k < j; k++)
                        {
                            if ((pts[j] - pts[k]).norm() < clamptol)
                                ok = false;
                        }
                        if (ok)
                            tokeep.push_back(j);
                    }

                    if (tokeep.size() > 2)
                    {
                        std::cerr << "warning: bad tetrahedron " << i << std::endl;
                    }
                    for (int pt1 = 0; pt1 < tokeep.size(); pt1++)
                    {
                        for (int pt2 = pt1 + 1; pt2 < tokeep.size(); pt2++)
                        {
                            Eigen::Vector3d ambpt[2];
                            ambpt[0].setZero();
                            ambpt[1].setZero();
                            for (int j = 0; j < 4; j++)
                            {
                                ambpt[0] += pts[tokeep[pt1]][j] * V.row(mesh.tetVertex(i, j));
                                ambpt[1] += pts[tokeep[pt2]][j] * V.row(mesh.tetVertex(i, j));
                            }
                            Segment s;
                            s.end1 = ambpt[0];
                            s.end2 = ambpt[1];
                            segs.push_back(s);
                        }
                    }
                }
            }
        }
    }

    std::cout << "BAD VERTS SIZE " << badverts.size() <<std::endl;

    int nsegs = segs.size();
    P.resize(2 * nsegs, 3);
    E.resize(nsegs, 2);
    for (int i = 0; i < nsegs; i++)
    {
        P.row(2 * i) = segs[i].end1;
        P.row(2 * i + 1) = segs[i].end2;
        E(i, 0) = 2 * i;
        E(i, 1) = 2 * i + 1;
    }

    int nsegs2 = segs2.size();
    P2.resize(2 * nsegs2, 3);
    E2.resize(nsegs2, 2);
    for (int i = 0; i < nsegs2; i++)
    {
        P2.row(2 * i) = segs2[i].end1;
        P2.row(2 * i + 1) = segs2[i].end2;
        E2(i, 0) = 2 * i;
        E2(i, 1) = 2 * i + 1;
    }

}