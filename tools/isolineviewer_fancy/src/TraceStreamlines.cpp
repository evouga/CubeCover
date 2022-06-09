#include "TraceStreamlines.h"
#include "TetMeshConnectivity.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "polyscope/point_cloud.h"

#include "FrameField.h"

#include <set>


static bool intersectPlane(const Eigen::Vector3d &n, const Eigen::Vector3d &p0, const Eigen::Vector3d &l0, const Eigen::Vector3d &l, double &t) 
{ 
    // assuming vectors are all normalized
    double denom = n.dot(l); 
    if (denom > 1e-6) { 
        Eigen::Vector3d p0l0 = p0 - l0; 
        t = p0l0.dot(n) / denom - .001; // in case tet is degenerate, don't quite hit the surface 
        return (t >= 0); 
    } 
 
    return false; 
}



static bool advanceStreamline( const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const CubeCover::FrameField& frames, 
    Streamline& s)
{
    int cur_tet_id = s.tetIds.back();
    Eigen::MatrixXd cur_tet_frame = frames.tetFrame(cur_tet_id);
    Eigen::Vector3d cur_point = s.points.back();
    Eigen::Vector3d prev_vec = s.directionVecs.back();

    double max_match = 0;
    Eigen::Vector3d first_frame_vec = cur_tet_frame.row(0);
    Eigen::Vector3d cur_vec = prev_vec.cross(first_frame_vec);
    for (int i = 0; i < cur_tet_frame.rows(); i++)
    {
        Eigen::Vector3d cur_dir = cur_tet_frame.row(i); 
        double cur_match = cur_dir.dot(prev_vec);
        if (std::abs(cur_match) > max_match)
        {
            max_match = std::abs(cur_match);
            cur_vec = cur_tet_frame.row(i);
            if (cur_match < 0)
                cur_vec = -cur_vec;
        }
    }
    cur_vec.normalize();


    // std::cout << "prev_vec: " << prev_vec  << std::endl;

    int cur_face_id = s.faceIds.back();
    int cur_local_fid = s.tetFaceIds.back();
    // int next_face_id = -1;

    // std::cout << "cur_vec: " << cur_vec;
    // std::cout << " cur_face_id: " << cur_face_id;
    // std::cout << " cur_local_fid: " << cur_local_fid ;
    // // std::cout << " cur_point: " << cur_point ;
    // std::cout << " cur_tet_id: " << cur_tet_id ;


    Eigen::MatrixXd face_normals;
    Eigen::MatrixXd face_basepoints;

    face_normals.resize(4,3);
    face_basepoints.resize(4,3);

    double min_intersect_time = 10000000;
    int next_tetface_id = -1;
    int next_face_id = -1;
    int next_tet_id = -1;
    Eigen::Vector3d intersect_point = cur_point;
    int next_face_orientation = -1;

    for (int i = 0; i < 4; i++)
    {
        if (i != cur_local_fid)
        {
            int fid = mesh.tetFace(cur_tet_id, i);
            Eigen::Vector3d v0 = V.row(mesh.faceVertex(fid, 0));
            Eigen::Vector3d v1 = V.row(mesh.faceVertex(fid, 1));
            Eigen::Vector3d v2 = V.row(mesh.faceVertex(fid, 2));

            Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
            n.normalize();
            if ( n.dot(cur_vec) < 0)
            {
                n = -n;
            } 
            face_normals.row(i) = n;
            face_basepoints.row(i) = v0;

            double intersect_time = 0;
            intersectPlane(n, v0, cur_point, cur_vec, intersect_time);
                // std::cout << "intersect_time: " << intersect_time  << std::endl;
            if (intersect_time < min_intersect_time && intersect_time > 0.000001)
            {
                min_intersect_time = std::abs(intersect_time);
                int cur_tetface_id = i;
                next_face_id = mesh.tetFace(cur_tet_id, cur_tetface_id);
                intersect_point = cur_point + cur_vec*intersect_time;
                next_face_orientation = mesh.tetFaceOrientation(cur_tet_id, i);
                next_face_orientation = (next_face_orientation + 1) % 2;
                next_tet_id = mesh.faceTet(next_face_id, next_face_orientation);
            }
        }
    }

    if ( next_face_id == -1)
    {
        std::cout << "Something wierd happened! Ray already at the boundary!" << std::endl;
        s.directionVecs.pop_back();
        s.points.pop_back();
        s.faceIds.pop_back();
        return false;
    }

    s.directionVecs.push_back(cur_vec);
    s.points.push_back(intersect_point);
    s.faceIds.push_back(next_face_id);

    // std::cout << "cur_local_fid: " << cur_local_fid;
    // std::cout << " next_tet_id: " << next_tet_id;
    // std::cout << " next_face_id: " << next_face_id ;
    // std::cout << " next_face_orientation: " << next_face_orientation ;
    // std::cout << " min_intersect_time: " << min_intersect_time ;

    if ( next_tet_id == -1)
    {
         // std::cout << "EXIT EXIT EXIT EXIT EXIT EXIT EXIT EXIT : " << next_face_orientation ;
        return false;
    }
    else
    {
        for (int i = 0; i < 4; i++)
        {
            int fid = mesh.tetFace(next_tet_id, i);
            if (fid == next_face_id)
            {
                next_tetface_id = i;
            }
        }
    }


    s.tetIds.push_back(next_tet_id);
    s.tetFaceIds.push_back(next_tetface_id);
    
    // std::cout << " next_tetface_id: " << next_tetface_id  << std::endl;


    return true;

}




void traceStreamlines(const Eigen::MatrixXd& V, const CubeCover::TetMeshConnectivity& mesh, const CubeCover::FrameField& frames,
    Eigen::VectorXi& init_tet_ids,
    int max_iter_per_trace,
    std::vector<Streamline>& traces)
{
    // constexpr double clamptol = 1e-6;

    // assert(values.rows() == mesh.nTets() * 4);
    // assert(values.cols() == 3);
    traces.clear();


    if (init_tet_ids.size() < 1)
    {
        std::vector<int> sampledTets;
        std::set<int> visitedTets;
        std::set<int>::iterator it;
        std::set<int>::iterator it2;
        std::pair<std::set<int>::iterator, bool> ret1;
        std::pair<std::set<int>::iterator, bool> ret2;

        int ntets = mesh.nTets();

        for (int i = 0; i < 100; i++)
        {
            it = visitedTets.find(i);
            if (it == visitedTets.end())
            {
                visitedTets.insert(i);
                sampledTets.push_back(i);
                // insert the neighbors.  
                for(int j = 0; j < 4; j++)
                {
                    int f1 = mesh.faceTet(mesh.tetFace(i, j), 0);
                    int f2 = mesh.faceTet(mesh.tetFace(i, j), 1);
                    ret1 = visitedTets.insert(f1);
                    ret2 = visitedTets.insert(f2);
                    it = ret1.first;
                    it2 = ret2.first;
                    int boundtet = -2;
                    if ( it == visitedTets.end())
                    {
                        if (it2 == visitedTets.end())
                        {
                            boundtet = -1;
                        }
                        else
                        {
                            boundtet = *it2;
                        }
                    }
                    else 
                    {
                        boundtet = *it;
                    }

                    if (boundtet > -1)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            f1 = mesh.faceTet(mesh.tetFace(boundtet, k), 0);
                            f2 = mesh.faceTet(mesh.tetFace(boundtet, k), 1);
                            ret1 = visitedTets.insert(f1);
                            ret2 = visitedTets.insert(f2);
                        }
                    }


                }

            }
        }
        int samplesize = sampledTets.size();
        init_tet_ids.resize(samplesize);
        for(int k = 0; k < samplesize; k++ )
        {
            init_tet_ids(k) = sampledTets.at(k);
        }

    }
    
    // std::vector<Streamline> traces2; 

    // std::vector<Segment> segs;
    // std::vector<Segment> segs2;1

    int counter = 0;

    int nstreamlines = init_tet_ids.rows();
    for (int i = 0; i < nstreamlines; i++)
    // for (int i = 20; i < 28; i++)
    {


        int cur_start_tet_id = init_tet_ids(i); 
        Eigen::MatrixXd cur_tet_frame = frames.tetFrame(cur_start_tet_id);
        int nvecs = cur_tet_frame.rows();
        
        // nvecs = 1; // hack
        for (int vec_id = 0; vec_id < nvecs; vec_id++)
        // for (int vec_id = 30; vec_id < 33; vec_id++)
        {
            Eigen::Vector3d cur_direction_vec = cur_tet_frame.row(vec_id);

            for (int orientation = 0; orientation < 2; orientation++)
            {
                if (orientation == 1)
                {
                    cur_direction_vec = -cur_direction_vec;
                }
                std::cout << std::endl << "processing streamline: " << counter << std::endl;
                counter++;


                Streamline t;
                t.tetIds.push_back(cur_start_tet_id);
                t.tetFaceIds.push_back(-1);

                int cur_face_id = -1;
                t.faceIds.push_back(cur_face_id);

                
                t.directionVecs.push_back(cur_direction_vec);


                Eigen::Vector4i cur_tet;
                cur_tet(0) = mesh.tetVertex(cur_start_tet_id, 0);
                cur_tet(1) = mesh.tetVertex(cur_start_tet_id, 1);
                cur_tet(2) = mesh.tetVertex(cur_start_tet_id, 2);
                cur_tet(3) = mesh.tetVertex(cur_start_tet_id, 3);

                Eigen::Vector3d cur_point;
                Eigen::Vector3d v0 = V.row(cur_tet(0));
                Eigen::Vector3d v1 = V.row(cur_tet(1));
                Eigen::Vector3d v2 = V.row(cur_tet(2));
                Eigen::Vector3d v3 = V.row(cur_tet(3));
                cur_point = ( v0 + v1 + v2 + v3 ) / 4.;
                t.points.push_back(cur_point);

                // Eigen::Vector3d cur_vec = cur_tet_frame.row(0);

                bool canContinue = true;

                for (int j = 0; j < max_iter_per_trace; j++)
                {
                    if (canContinue)
                    {
                        canContinue = advanceStreamline(V, mesh, frames, t);
                    }
                    else
                    {
                        break;
                    }
                }
                if (t.tetIds.size() > 1)
                {
                    traces.push_back(t);
                }
                
            }
        }
    }

}

    //     // Eigen::MatrixXd cur_face;


    //     double minvals[3];
    //     double maxvals[3];
    //     for (int j = 0; j < 3; j++)
    //     {
    //         minvals[j] = std::numeric_limits<double>::infinity();
    //         maxvals[j] = -std::numeric_limits<double>::infinity();
    //     }
    //     for (int j = 0; j < 4; j++)
    //     {
    //         for (int k = 0; k < 3; k++)
    //         {
    //             minvals[k] = std::min(minvals[k], values(4 * i + j, k));
    //             maxvals[k] = std::max(maxvals[k], values(4 * i + j, k));
    //         }
    //     }

    //     int tosamplelower[3];
    //     int tosampleupper[3];

    //     for (int j = 0; j < 3; j++)
    //     {
    //         tosamplelower[j] = std::ceil(minvals[j]);
    //         tosampleupper[j] = std::floor(maxvals[j]);
    //     }

    //     // add segments for each direction
    //     for (int dir = 0; dir < 3; dir++)
    //     {            
    //         int sample1 = (dir + 1) % 3;
    //         int sample2 = (dir + 2) % 3;
    //         for (int val1 = tosamplelower[sample1]; val1 <= tosampleupper[sample1]; val1++)
    //         {
    //             for (int val2 = tosamplelower[sample2]; val2 <= tosampleupper[sample2]; val2++)
    //             {
    //                 std::vector<Eigen::Vector4d> pts;
    //                 int nvertpts = 0;
    //                 int nedgepts = 0;
    //                 int nfacepts = 0;

    //                 // vertex case
    //                 for (int vert = 0; vert < 4; vert++)
    //                 {
    //                     double s1val = values(4 * i + vert, sample1);
    //                     double s2val = values(4 * i + vert, sample2);
    //                     if (s1val == double(val1) && s2val == double(val2))
    //                     {
    //                         Eigen::Vector4d pt(0, 0, 0, 0);
    //                         pt[vert] = 1.0;
    //                         pts.push_back(pt);
    //                         nvertpts++;
    //                     }
    //                 }

    //                 // edge case
    //                 for (int vert1 = 0; vert1 < 4; vert1++)
    //                 {
    //                     for (int vert2 = vert1 + 1; vert2 < 4; vert2++)
    //                     {
    //                         double s11val = values(4 * i + vert1, sample1);
    //                         double s12val = values(4 * i + vert2, sample1);
    //                         double s21val = values(4 * i + vert1, sample2);
    //                         double s22val = values(4 * i + vert2, sample2);
    //                         if (s11val == double(val1) && s12val == double(val1))
    //                         {
    //                             // whole edge intersect the isosurface
    //                             if (s21val != s22val)
    //                             {
    //                                 // equality is handled in the vertex case
    //                                 double bary = (double(val2) - s21val) / (s22val - s21val);
    //                                 if (bary > 0.0 && bary < 1.0)
    //                                 {                                        
    //                                     Eigen::Vector4d pt(0, 0, 0, 0);
    //                                     pt[vert1] = 1.0 - bary;
    //                                     pt[vert2] = bary;
    //                                     pts.push_back(pt);
    //                                     nedgepts++;
    //                                 }
    //                             }
    //                             else
    //                             {
    //                                 double bary1 = (double(val1) - s11val) / (s12val - s11val);
    //                                 double bary2 = (double(val2) - s21val) / (s22val - s21val);
    //                                 if (bary1 > 0.0 && bary1 < 1.0 && bary1 == bary2)
    //                                 {
    //                                     Eigen::Vector4d pt(0, 0, 0, 0);
    //                                     pt[vert1] = 1.0 - bary1;
    //                                     pt[vert2] = bary1;
    //                                     pts.push_back(pt);
    //                                     nedgepts++;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }


    //                 // face case
    //                 for (int vert1 = 0; vert1 < 4; vert1++)
    //                 {
    //                     for (int vert2 = vert1 + 1; vert2 < 4; vert2++)
    //                     {
    //                         for (int vert3 = vert2 + 1; vert3 < 4; vert3++)
    //                         {
    //                             double s11val = values(4 * i + vert1, sample1);
    //                             double s12val = values(4 * i + vert2, sample1);
    //                             double s13val = values(4 * i + vert3, sample1);
    //                             double s21val = values(4 * i + vert1, sample2);
    //                             double s22val = values(4 * i + vert2, sample2);
    //                             double s23val = values(4 * i + vert3, sample2);
                                
    //                             Eigen::Matrix2d M;
    //                             M << s12val - s11val, s13val - s11val,
    //                                 s22val - s21val, s23val - s21val;

    //                             if (M.determinant() == 0)
    //                             {
    //                                 // std::cout << "zero volume" <<std::endl;
    //                                 // Segment s;
    //                                 // s.end1 = V.row( mesh.tetVertex(i, vert1) );
    //                                 // s.end2 = V.row( mesh.tetVertex(i, vert2) );;
    //                                 // segs2.push_back(s);
    //                                 // s.end1 = V.row( mesh.tetVertex(i, vert2) );
    //                                 // s.end2 = V.row( mesh.tetVertex(i, vert3) );;
    //                                 // segs2.push_back(s);
    //                                 // s.end1 = V.row( mesh.tetVertex(i, vert1) );
    //                                 // s.end2 = V.row( mesh.tetVertex(i, vert3) );;
    //                                 // segs2.push_back(s);
    //                                 continue;
    //                             }

    //                             Eigen::Vector2d rhs(double(val1) - s11val, double(val2) - s21val);
    //                             Eigen::Vector2d barys = M.inverse() * rhs;
    //                             if (barys[0] >= 0 && barys[1] >= 0 && barys[0] + barys[1] <= 1.0)
    //                             {
    //                                 Eigen::Vector4d pt(0, 0, 0, 0);
    //                                 pt[vert1] = 1.0 - barys[0] - barys[1];
    //                                 pt[vert2] = barys[0];
    //                                 pt[vert3] = barys[1];
    //                                 pts.push_back(pt);
    //                                 nfacepts++;
    //                             }                                
    //                         }
    //                     }
    //                 }

    //                 // clamp together identical pts
    //                 std::vector<int> tokeep;
    //                 int ncrossings = pts.size();
    //                 for (int j = 0; j < ncrossings; j++)
    //                 {
    //                     bool ok = true;
    //                     for (int k = 0; k < j; k++)
    //                     {
    //                         if ((pts[j] - pts[k]).norm() < clamptol)
    //                             ok = false;
    //                     }
    //                     if (ok)
    //                         tokeep.push_back(j);
    //                 }

    //                 if (tokeep.size() > 2)
    //                 {
    //                     std::cerr << "warning: bad tetrahedron " << i << std::endl;
    //                 }
    //                 for (int pt1 = 0; pt1 < tokeep.size(); pt1++)
    //                 {
    //                     for (int pt2 = pt1 + 1; pt2 < tokeep.size(); pt2++)
    //                     {
    //                         Eigen::Vector3d ambpt[2];
    //                         ambpt[0].setZero();
    //                         ambpt[1].setZero();
    //                         for (int j = 0; j < 4; j++)
    //                         {
    //                             ambpt[0] += pts[tokeep[pt1]][j] * V.row(mesh.tetVertex(i, j));
    //                             ambpt[1] += pts[tokeep[pt2]][j] * V.row(mesh.tetVertex(i, j));
    //                         }
    //                         Segment s;
    //                         s.end1 = ambpt[0];
    //                         s.end2 = ambpt[1];
    //                         segs.push_back(s);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // int nsegs = segs.size();
    // P.resize(2 * nsegs, 3);
    // E.resize(nsegs, 2);
    // for (int i = 0; i < nsegs; i++)
    // {
    //     P.row(2 * i) = segs[i].end1;
    //     P.row(2 * i + 1) = segs[i].end2;
    //     E(i, 0) = 2 * i;
    //     E(i, 1) = 2 * i + 1;
    // }

    // int nsegs2 = segs2.size();
    // P2.resize(2 * nsegs2, 3);
    // E2.resize(nsegs2, 2);
    // for (int i = 0; i < nsegs2; i++)
    // {
    //     P2.row(2 * i) = segs2[i].end1;
    //     P2.row(2 * i + 1) = segs2[i].end2;
    //     E2(i, 0) = 2 * i;
    //     E2(i, 1) = 2 * i + 1;
    // }

// }