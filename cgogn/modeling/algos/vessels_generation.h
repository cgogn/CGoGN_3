/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#ifndef CGOGN_GEOMETRY_ALGOS_VESSELS_GENERATION_H_
#define CGOGN_GEOMETRY_ALGOS_VESSELS_GENERATION_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/functions/intersection.h>


namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;
using Eigen::Matrix3d;


using GDart = Graph::Vertex1;
using GVertex = Graph::Vertex;
using GEdge = Graph::Edge;

using M2Vertex = CMap2::Vertex;
using M2Edge = CMap2::Edge;
using M2Face = CMap2::Face;
using M2Volume = CMap2::Volume;

using GBranch = std::pair<Dart, Dart>;

struct GraphData
{
    std::vector<GBranch> branches;
    std::vector<GVertex> intersections;
};


struct GraphAttributes
{
    std::shared_ptr<Graph::Attribute<geometry::Vec3>> V_pos;
    std::shared_ptr<Graph::Attribute<geometry::Scalar>> V_rad;
    std::shared_ptr<Graph::Attribute<Dart>> m2_CC;
    std::shared_ptr<Graph::Attribute<Dart>> m2_interface;
    std::shared_ptr<Graph::Attribute<Matrix3d>> frames;
    
};

struct M2Attributes
{
    std::shared_ptr<CMap2::Attribute<geometry::Vec3>> V_pos;
    std::shared_ptr<CMap2::Attribute<geometry::Vec3>> center;
    std::shared_ptr<CMap2::Attribute<geometry::Vec3>> f_point;
    std::shared_ptr<CMap2::Attribute<Dart>> f_branch;
};

bool get_graph_data(Graph& g, GraphData& gData)
{
    std::vector<GVertex> extremities;

    foreach_cell(g, [&](GVertex gv) -> bool 
    {
        switch(nb_darts_of_orbit(g, gv)){
            case 1:
                extremities.push_back(gv);
                break;
            case 2:
                break;
            default:
                extremities.push_back(gv);
                gData.intersections.push_back(gv);
                break;
        }
        return true;
    });

    CellMarker<Graph, GDart> cmg(g);
    for(GVertex gv : extremities)
    {
        g.foreach_dart_of_orbit(gv, [&](Dart d0)-> bool 
        {
            if(!cmg.is_marked(GDart(d0)))
            {
                cmg.mark(GDart(d0));
                Dart d1 = g.alpha0(d0);
                cmg.mark(GDart(d1));
                while(nb_darts_of_orbit(g, GVertex(d1)) == 2)
                {
                    d1 = g.alpha1(d1);
                    cmg.mark(GDart(d1));
                    d1 = g.alpha0(d1);
                    cmg.mark(GDart(d1));
                }
                gData.branches.push_back(std::pair(d0, d1));
            }
            return true;
        });
    }

    return true;
}

bool get_graph_attributes(Graph& g, GraphAttributes& gAttribs)
{
    gAttribs.V_pos = cgogn::get_attribute<geometry::Vec3, Graph::Vertex>(g, "position");
    if(!gAttribs.V_pos)
    {
        std::cout << "The graph has no position atttribute" << std::endl;
        return false;
    }

    gAttribs.V_rad = cgogn::get_attribute<geometry::Scalar, GVertex>(g, "radius");
    if(!gAttribs.V_rad)
    {
        std::cout << "The                 create_intersection_frames(g, gAttribs, m2, m2Attribs, gv);graph has no radius atttribute" << std::endl;
        return false;
    }

    gAttribs.m2_CC = cgogn::add_attribute<Dart, GVertex>(g, "m2_CC");
    if(!gAttribs.m2_CC)
    {
        std::cout << "Failed to add m2_CC attribute to graph" << std::endl;
        return false;
    }

    gAttribs.m2_interface = cgogn::add_attribute<Dart, GDart>(g, "m2_interface");
    if(!gAttribs.m2_interface)
    {
        std::cout << "Failed to add m2_CC attribute to graph" << std::endl;
        return false;
    }

    gAttribs.frames = cgogn::add_attribute<Matrix3d, GDart>(g, "frames");
    if(!gAttribs.frames)
    {
        std::cout << "Failed to add Frames attribute to graph" << std::endl;
        return false;
    }

    return true;
}

// bool clean_graph(Graph& g)
// {
//     return true;
// }

bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs)
{
    m2Attribs.V_pos = cgogn::add_attribute<Vec3, M2Vertex>(m2, "position");
    if(!m2Attribs.V_pos)
    {
        std::cout << "Failed to add position attribute to map2" << std::endl;
        return false;
    }

    m2Attribs.center = cgogn::add_attribute<Vec3, M2Volume>(m2, "center");
    if(!m2Attribs.center)
    {
        std::cout << "Failed to add center attribute to map2" << std::endl;
        return false;
    }

    m2Attribs.f_point = cgogn::add_attribute<Vec3, M2Face>(m2, "f_point");
    if(!m2Attribs.f_point)
    {
        std::cout << "Failed to add f_point attribute to map2" << std::endl;
        return false;
    }

    m2Attribs.f_branch = cgogn::add_attribute<Dart, M2Face>(m2, "f_branch");
    if(!m2Attribs.f_branch)
    {
        std::cout << "Failed to add f_branch attribute to map2" << std::endl;
        return false;
    }

    return true;
}

Vec3 project_on_sphere(Vec3 P, Scalar R, Vec3 C)
{
    return C + (P - C).normalized() * R;
}

bool build_interface_1(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, GVertex gv)
{
    M2Face f = add_face(m2, 4, false);

    value<Dart>(g, gAttribs.m2_CC, gv) = f.dart;
    value<Dart>(g, gAttribs.m2_interface, GDart(gv.dart)) = f.dart;

    return true;
}

bool build_interface_2(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, GVertex gv)
{
    Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
    Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;

    m2.phi2_sew(d0, d1);
    m2.phi2_sew(m2.phi1(d0), m2.phi_1(d1));
    m2.phi2_sew(m2.phi<11>(d0), m2.phi<11>(d1));
    m2.phi2_sew(m2.phi_1(d0), m2.phi1(d1));

    value<Dart>(g, gAttribs.m2_CC, GDart(gv.dart)) = d0;
    value<Dart>(g, gAttribs.m2_interface, GDart(gv.dart)) = d0;
    value<Dart>(g, gAttribs.m2_interface, GDart(g.alpha1(gv.dart))) = d1;

    return true;
}

bool build_interface_core(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, GVertex gv)
{

    M2Face f0 = add_face(static_cast<CMap1&>(m2), 4, false);
    M2Face f1 = add_face(static_cast<CMap1&>(m2), 4, false);
    M2Face f2 = add_face(static_cast<CMap1&>(m2), 4, false);

    Dart d0 = f0.dart;
    Dart d1 = f1.dart;
    Dart d2 = f2.dart;

    m2.phi2_sew(d0, m2.phi1(d1));
    m2.phi2_sew(m2.phi1(d0), d2);
    m2.phi2_sew(m2.phi<11>(d0), m2.phi_1(d2));
    m2.phi2_sew(m2.phi_1(d0), m2.phi<11>(d1));
    m2.phi2_sew(d1, m2.phi1(d2));
    m2.phi2_sew(m2.phi_1(d1), m2.phi<11>(d2));

    value<Dart>(g, gAttribs.m2_CC, gv) = d0;
    return true;
}

bool build_connection_interfaces(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2)
{
    foreach_cell(g, [&](GVertex gv) -> bool
    {
        switch(nb_darts_of_orbit(g, gv))
        {
            case 1:
                build_interface_1(g, gAttribs, m2, gv);
                break;
            case 2:
                build_interface_2(g, gAttribs, m2, gv);
                break;
            case 3:
            default:
                build_interface_core(g, gAttribs, m2, gv);
                break;
        }
        return true;
    });
    return true;
}

bool complete_intersection_3(const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs, GVertex gv)
{
    
    Vec3 center = value<Vec3>(g, gAttribs.V_pos, gv);
    Scalar radius = value<Scalar>(g, gAttribs.V_rad, gv);
    Dart cc = value<Dart>(g, gAttribs.m2_CC, gv); 

    std::vector<Dart> Fd = {cc, m2.phi<2111>(cc), m2.phi<12>(cc)};

    std::vector<Vec3> Ppos;
    Ppos.reserve(3);
    std::vector<Dart> Pdart;
    Pdart.reserve(3);

    g.foreach_dart_of_orbit(gv, [&](Dart d) -> bool 
    {
        Vec3 p = value<Vec3>(g, gAttribs.V_pos, GVertex(g.alpha0(d)));
        p = project_on_sphere(p, radius, center);
        Ppos.push_back(p);
        Pdart.push_back(d);
        return true;
    });

    Vec3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
    std::vector<Vec3> Q = {center + V * radius, center - V * radius};
    std::vector<Vec3> M = {center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
        center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
        center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};
    
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Fd[0])) = M[0];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Fd[1])) = M[1];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Fd[2])) = M[2];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi1(Fd[0]))) = Q[0];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi_1(Fd[0]))) = Q[1];

    value<Dart>(g, gAttribs.m2_interface, GDart(Pdart[0])) = Fd[0];
    value<Dart>(g, gAttribs.m2_interface, GDart(Pdart[1])) = Fd[1];
    value<Dart>(g, gAttribs.m2_interface, GDart(Pdart[2])) = Fd[2];

    return true;
}

Dart convex_quad(CMap2& m2, M2Attributes& m2Attribs, Dart f){
        Dart res;
        Vec3 A = value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(f));
        Vec3 B = value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi1(f))); 
        Vec3 C = value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi<11>(f)));
        Vec3 D = value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi_1(f))); 
        Vec3 AB = B - A;
        Vec3 AC = C - A;
        Vec3 AD = D - A;
        Vec3 N0 = AC.cross(AD);
        Vec3 N1 = AB.cross(AC);
        Vec3 N = N0.cross(N1).normalized();

        if(N.dot(AC) < 0){
            res = m2.phi1(f);
        }
        else{
            res = f;
        }
        return res;
    }

Vec3 mean_dir(Vec3 center, Scalar radius, Vec3 point, std::vector<Vec3> points){
    using Quat = Eigen::Quaterniond;
    uint32 valence = points.size();

    std::vector<Vec3> directions;
    for(Vec3 p : points){
        directions.push_back((p - center).normalized());
    }
    Vec3 avg_dir = (point - center).normalized();

    std::vector<Quat> rotations;
    rotations.reserve(valence);

//        for(uint i = 0; i < valence; i++){
        for(Vec3 dir : directions){
            Quat q = Quat::FromTwoVectors(avg_dir, dir);
            q.normalize();
            rotations.push_back(q);
        }

        Eigen::MatrixXd m(4, valence);
        for(uint32 j = 0; j < valence; ++j){
            const Quat& q = rotations[j];
            m.col(j) = Eigen::Vector4d(q.w(), q.x(), q.y(), q.z());
        }

        Eigen::MatrixXd mm = m * m.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mm);
        Eigen::Vector4d r = es.eigenvectors().col(3);

        Quat mean_rot(r[0], r[1], r[2], r[3]);
        mean_rot.normalize();
        avg_dir = mean_rot._transformVector(avg_dir);
        rotations.clear();
//        }

    return avg_dir * radius + center;
}

bool complete_intersection_n(const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs, GVertex gv)
{
    using cgogn::geometry::intersection_ray_triangle;

    Vec3 center = value<Vec3>(g, gAttribs.V_pos, gv);
    Scalar radius = value<Scalar>(g, gAttribs.V_rad, gv);
    Dart cc = value<Dart>(g, gAttribs.m2_CC, gv); 

    std::vector<Dart> Fd = {cc, m2.phi<2111>(cc), m2.phi<12>(cc)};
    std::vector<Vec3> Ppos;
    Ppos.reserve(3);
    std::vector<Dart> Pdart;
    Pdart.reserve(3);
    g.foreach_dart_of_orbit(gv, [&](Dart d) -> bool 
    {
        Vec3 p = value<Vec3>(g, gAttribs.V_pos, GVertex(g.alpha0(d)));
        p = project_on_sphere(p, radius, center);
        Ppos.push_back(p);
        Pdart.push_back(d);
        return true;
    });

    Vec3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
    std::vector<Vec3> Q = {center + V * radius, center - V * radius};
    std::vector<Vec3> M = {center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
        center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
        center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};
    
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Fd[0])) = M[0];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Fd[1])) = M[1];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Fd[2])) = M[2];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi1(Fd[0]))) = Q[0];
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi_1(Fd[0]))) = Q[1];

    value<Vec3>(m2, m2Attribs.f_point, M2Face(Fd[0])) = Ppos[0];
    value<Vec3>(m2, m2Attribs.f_point, M2Face(Fd[1])) = Ppos[1];
    value<Vec3>(m2, m2Attribs.f_point, M2Face(Fd[2])) = Ppos[2];

    value<Dart>(m2, m2Attribs.f_branch, M2Face(Fd[0])) = Pdart[0];
    value<Dart>(m2, m2Attribs.f_branch, M2Face(Fd[1])) = Pdart[1];
    value<Dart>(m2, m2Attribs.f_branch, M2Face(Fd[2])) = Pdart[2];

    for(uint i = 3; i < Ppos.size(); i++)
    {
        Vec3 P0 = Ppos[i];
        Dart F0 = Pdart[i];
        M2Face face2cut;

        std::vector<Vec3> Quadp;
        std::vector<Dart> Quadd;

        bool face_found = false;

        foreach_cell(m2, [&](M2Face f) -> bool 
        {
            Dart fd = convex_quad(m2, m2Attribs, f.dart);
            Quadd = {fd, m2.phi1(fd),
                        m2.phi<11>(fd),
                        m2.phi_1(fd)};
            Quadp = {value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Quadd[0])),
                    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Quadd[1])),
                    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Quadd[2])),
                    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(Quadd[3]))};
            if((face_found = (intersection_ray_triangle(center, P0 - center, Quadp[0], Quadp[1], Quadp[2])
                                || intersection_ray_triangle(center, P0 - center, Quadp[0], Quadp[2], Quadp[3])))){
                face2cut = f;
                return false;
            }
            return true;
        });
        if(!face_found)
        {
            std::cout << "error complete_intersection_n: no face found" << std::endl;
            return false;
        }

        Vec3 P1 = value<Vec3>(m2, m2Attribs.f_point, M2Face(face2cut));
        Dart F1 = value<Dart>(m2, m2Attribs.f_branch, M2Face(face2cut));
        
        Dart cut0, cut1;
        Vec3 AC = (Quadp[2] - Quadp[0]).normalized();
        Vec3 BD = (Quadp[3] - Quadp[1]).normalized();
        Vec3 P0P1 = (P1 - P0).normalized();
        if(abs(AC.dot(P0P1)) < abs(BD.dot(P0P1))){
            cut0 = Quadd[0]; cut1 = Quadd[2];
        } else {
            cut0 = Quadd[1]; cut1 = Quadd[3];
        }
        cut_face(m2, M2Vertex(cut0), M2Vertex(cut1));
        M2Vertex v = cut_edge(m2, M2Edge(m2.phi_1(cut0)));
        value<Vec3>(m2, m2Attribs.V_pos, v) = project_on_sphere((P0 + P1) * Scalar(0.5), radius, center);

        Dart newFace, oldFace;
        Vec3 out0 =  value<Vec3>(m2, m2Attribs.V_pos, v) - center;
        // Vec3 out1 = ((P0 - value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(cut0))).normalized().cross((P1 - m2_attribs.position[cut0]).normalized()));

        Vec3 out1 = ((P0 - value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(cut0))).normalized().cross((P1 - value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(cut0))).normalized()));
        if(out1.dot(out0) >= 0)
        {
            newFace = cut0; oldFace = cut1;
        }
        else
        {
            newFace = cut1; oldFace = cut0;
        }

        value<Vec3>(m2, m2Attribs.f_point, M2Face(newFace)) = P0;
        value<Vec3>(m2, m2Attribs.f_point, M2Face(oldFace)) = P1;
        value<Dart>(m2, m2Attribs.f_branch, M2Face(newFace)) = F0;
        value<Dart>(m2, m2Attribs.f_branch, M2Face(oldFace)) = F1;

        foreach_incident_vertex(m2, M2Face(newFace), [&](M2Vertex m2v) -> bool
        {
            std::vector<Vec3> points;
            foreach_incident_face(m2, m2v, [&](M2Face m2f) -> bool
            {
                points.push_back(value<Vec3>(m2, m2Attribs.f_point, m2f));
                return true;
            });
            value<Vec3>(m2, m2Attribs.V_pos, m2v) = mean_dir(center, radius, value<Vec3>(m2, m2Attribs.V_pos, m2v), points);
            return true;
        });
    }

    foreach_incident_face(m2, M2Volume(cc), [&](M2Face m2f) -> bool
    {
        Dart f = value<Dart>(m2, m2Attribs.f_branch, m2f);
        value<Dart>(g, gAttribs.m2_interface, GDart(f)) = m2f.dart;
        return true;
    });

    return true;
}

bool complete_intersections(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs)
{
    foreach_cell(g, [&](GVertex gv) -> bool
    {
        switch(nb_darts_of_orbit(g, gv))
        {
            case 1:
                break;
            case 2:
                break;
            case 3:
                complete_intersection_3(g, gAttribs, m2, m2Attribs, gv);
                break;
            default:
                complete_intersection_n(g, gAttribs, m2, m2Attribs, gv);
                break;
        }
        return true;
    });

    return true;
}

bool create_intersection_frames(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs, GVertex gv)
{
    Vec3 center = value<Vec3>(g, gAttribs.V_pos, gv);
    g.foreach_dart_of_orbit(gv, [&](Dart d) -> bool
    {
        std::cout << "dart :" << index_of(g, GDart(d)) << std::endl;
        Dart d0 = value<Dart>(g, gAttribs.m2_interface, GDart(d));
        Dart d1 = m2.phi<11>(d0);
        std::cout << "M2_interface :" << d0 << std::endl;

        Vec3 R, S, T, diag, temp;
        T = (value<Vec3>(g, gAttribs.V_pos, GVertex(g.alpha0(d))) - center).normalized();
        diag =(value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(d1)) -value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(d0))).normalized();
        R = diag.cross(T).normalized();
        S = T.cross(R).normalized();

        value<Matrix3d>(g, gAttribs.frames, GDart(d)).col(0) = R;
        value<Matrix3d>(g, gAttribs.frames, GDart(d)).col(1) = S;
        value<Matrix3d>(g, gAttribs.frames, GDart(d)).col(2) = T;

        return true;
    });
    return true;
}

bool create_intersections_frames(
    const Graph& g, GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs)
{
    foreach_cell(g, [&](GVertex gv) -> bool
    {
        switch(nb_darts_of_orbit(g, gv))
        {
            case 1:
                break;
            case 2:
                break;
            default:
                create_intersection_frames(g, gAttribs, m2, m2Attribs, gv);
                break;
        }
        return true;
    });

    return true;
}

Dart shift_interface(CMap2& m2, Dart m2d, uint32 nb_shifts){
        Dart d = m2d;
        for(uint i = 0; i < nb_shifts; ++i)
            d = m2.phi1(d);
        return d;
    }

Matrix3d shift_frame(Matrix3d frame, uint nb_shifts){
    Matrix3d f = frame;
    for(uint i = 0; i < nb_shifts; ++i){
        Vec3 R, S, T;
        R = f.col(1);
        S = -f.col(0);
        T = f.col(2);
        f.col(0) = R;
        f.col(1) = S;
        f.col(2) = T;
    }
    return f;
}

bool propagate_frame_n_n(const Graph& g, GraphAttributes& gAttribs, CMap2& m2, Dart d)
{
    Matrix3d U0 = value<Matrix3d>(g, gAttribs.frames, GDart(d));
    Dart d0 = d;
    Dart d1 = g.alpha0(d0);
    uint32 valence = nb_darts_of_orbit(g, GVertex(d1));
    Matrix3d U, U_;
    U = U0;

    uint nb_e = 0;

    Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
    Scalar c1, c2;

    while(valence == 2){
        nb_e++;

        Ri = U.col(0);
        Ti = U.col(2);
        v1 = value<Vec3>(g, gAttribs.V_pos, GVertex(d1)) - value<Vec3>(g, gAttribs.V_pos, GVertex(d0));
        c1 = v1.dot(v1);
        RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
        TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;

        v3 = (value<Vec3>(g, gAttribs.V_pos, GVertex(g.alpha0(g.alpha1(d1)))) - value<Vec3>(g, gAttribs.V_pos, GVertex(d1)));
        Ti1 = (v1.normalized()+v3.normalized()).normalized();
        v2 = Ti1 - TiL;
        c2 = v2.dot(v2);
        Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
        Si1 = Ti1.cross(Ri1);
        U.col(0) = Ri1;
        U.col(1) = Si1;
        U.col(2) = Ti1;

        d0 = g.alpha1(d1);
        d1 = g.alpha0(d0);
        value<Matrix3d>(g, gAttribs.frames, GDart(d0)) = U;
        valence = nb_darts_of_orbit(g, GVertex(d1));
    }

    Ri = U.col(0);
    Ti = U.col(2);
    v1 = value<Vec3>(g, gAttribs.V_pos, GVertex(d1)) - value<Vec3>(g, gAttribs.V_pos, GVertex(d0));
    c1 = v1.dot(v1);
    RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
    TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;
    Ti1 = v1.normalized();
    v2 = Ti1 - TiL;
    c2 = v2.dot(v2);
    Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
    Si1 = Ti1.cross(Ri1);
    U.col(0) = Ri1;
    U.col(1) = Si1;
    U.col(2) = Ti1;

    U_.col(0) = U.col(0);
    U_.col(1) = -U.col(1);
    U_.col(2) = -U.col(2);

    Vec3 X = (U_.col(0) + U_.col(1)).normalized();
    Matrix3d UE = value<Matrix3d>(g, gAttribs.frames, GDart(d1));
    Vec3 RE = UE.col(0), SE = UE.col(1);
    bool A = (RE.dot(X) >= 0);
    bool B = (SE.dot(X) >= 0);
    uint nb_shifts = 0;
    if(!A && B) nb_shifts = 1;
    else if(!A && !B) nb_shifts = 2;
    else if(A && !B) nb_shifts = 3;
    if(nb_shifts){
        UE = shift_frame(UE, nb_shifts);
        value<Matrix3d>(g, gAttribs.frames, GDart(d1)) = UE;
        value<Dart>(g, gAttribs.m2_interface, GDart(d1)) = shift_interface(m2, value<Dart>(g, gAttribs.m2_interface, GDart(d1)), nb_shifts);
    }

    if(nb_e){
        Scalar cos = UE.col(0).dot(U_.col(0));
        Scalar angle = cos > 1 ? std::acos(1) : std::acos(cos);
        Scalar angle_step = angle / nb_e;

        Dart d0 = d;
        Dart d1 = g.alpha0(d0);
        uint valence = nb_darts_of_orbit(g, GVertex(d1));
        uint step = 0;

        while(valence == 2) {
            step++;
            d0 = g.alpha1(d1);
            U = value<Matrix3d>(g, gAttribs.frames, GDart(d0));
            Eigen::AngleAxisd rot (angle_step * step, U.col(2));
            U.col(0) = rot * U.col(0);
            U.col(1) = U.col(2).cross(U.col(0));

            U_.col(0) = U.col(0);
            U_.col(1) = -U.col(1);
            U_.col(2) = -U.col(2);

            value<Matrix3d>(g, gAttribs.frames, GDart(d1)) = U_;
            value<Matrix3d>(g, gAttribs.frames, GDart(d0)) = U;

            d1 = g.alpha0(d0);
            valence = nb_darts_of_orbit(g, GVertex(d1));
        }
    }
    return true;
}

bool propagate_frame_n_1(const Graph& g, GraphAttributes& gAttribs, Dart d)
{
    Matrix3d U0 = value<Matrix3d>(g, gAttribs.frames, GDart(d));
    Dart d0 = d;
    Dart d1 = g.alpha0(d0);
    uint32 valence = nb_darts_of_orbit(g, GVertex(d1));
    Matrix3d U, U_;
    U = U0;

    while(valence == 2)
    {
        Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
        Scalar c1, c2;

        Ri = U.col(0);
        Ti = U.col(2);
        v1 = value<Vec3>(g, gAttribs.V_pos, GVertex(d1)) - value<Vec3>(g, gAttribs.V_pos, GVertex(d0));
        c1 = v1.dot(v1);
        RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
        TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;
        
        v3 = (value<Vec3>(g, gAttribs.V_pos, GVertex(g.alpha0(g.alpha1(d1)))) - value<Vec3>(g, gAttribs.V_pos, GVertex(d1)));
        Ti1 = (v1.normalized()+v3.normalized()).normalized();
        v2 = Ti1 - TiL;
        c2 = v2.dot(v2);
        Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
        Si1 = Ti1.cross(Ri1);

        U.col(0) = Ri1;
        U.col(1) = Si1;
        U.col(2) = Ti1;

        U_.col(0) = U.col(0);
        U_.col(1) = -U.col(1);
        U_.col(2) = -U.col(2);

        value<Matrix3d>(g, gAttribs.frames, GDart(d1)) = U_;
        d0 = g.alpha1(d1);
        d1 = g.alpha0(d0);
        value<Matrix3d>(g, gAttribs.frames, GDart(d0)) = U;
        valence = nb_darts_of_orbit(g, GVertex(d1));
    }
    if(valence == 1)
    {
        Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
        Scalar c1, c2;

        Ri = U.col(0);
        Ti = U.col(2);
        v1 = value<Vec3>(g, gAttribs.V_pos, GVertex(d1)) - value<Vec3>(g, gAttribs.V_pos, GVertex(d0));
        c1 = v1.dot(v1);
        RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
        TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;
        Ti1 = v1.normalized();
        v2 = Ti1 - TiL;
        c2 = v2.dot(v2);
        Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
        Si1 = Ti1.cross(Ri1);

        U.col(0) = Ri1;
        U.col(1) = Si1;
        U.col(2) = Ti1;

        U_.col(0) = U.col(0);
        U_.col(1) = -U.col(1);
        U_.col(2) = -U.col(2);

        value<Matrix3d>(g, gAttribs.frames, GDart(d1)) = U_;
    }
}

bool propagate_frames(const Graph& g, GraphAttributes& gAttribs, const GraphData& gData, CMap2& m2)
{
    for(auto branch : gData.branches)
    {
        if(nb_darts_of_orbit(g, GVertex(branch.first)) > 1)
        {
            if(nb_darts_of_orbit(g, GVertex(branch.second)) == 1)
                propagate_frame_n_1(g, gAttribs, branch.first);
            else
                propagate_frame_n_n(g, gAttribs, m2, branch.first);
        }
        else
        {
            propagate_frame_n_1(g, gAttribs, branch.second);
        }
    }
    return true;
}

bool complete_interface_1(const Graph& g, const GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs, GVertex gv)
{
    // Dart cc = value<Dart>(g, gAttribs.m2_CC, gv);
    Matrix3d frame = value<Matrix3d>(g, gAttribs.frames, GDart(gv.dart));
    Vec3 center = value<Vec3>(g, gAttribs.V_pos, gv);
    Scalar radius = value<Scalar>(g, gAttribs.V_rad, gv);
    Dart f0 = value<Dart>(g, gAttribs.m2_interface, GDart(gv.dart));

    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(f0)) = center - frame.col(1) * radius;
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi1(f0))) = center + frame.col(0) * radius;
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi<11>(f0))) = center + frame.col(1) * radius;
    value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(m2.phi_1(f0))) = center - frame.col(0) * radius;

    return true;
}

bool complete_interface_2(const Graph& g, const GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs, GVertex gv)
{
    Dart cc = value<Dart>(g, gAttribs.m2_CC, gv);
    m2.foreach_dart_of_orbit(M2Face(cc), [&](Dart d) -> bool 
    {
        value<Vec3>(m2, m2Attribs.V_pos, M2Vertex(d)) = value<Vec3>(g, gAttribs.V_pos, gv);
        return true;
    });

    return true;
}

bool complete_interfaces(const Graph& g, const GraphAttributes& gAttribs, 
    CMap2& m2, M2Attributes& m2Attribs)
{
    foreach_cell(g, [&](GVertex gv) -> bool
    {
        switch(nb_darts_of_orbit(g, gv))
        {
            case 1:
            case 2:
                complete_interface_1(g, gAttribs, m2, m2Attribs, gv);
                break;
            default:
                break;
        }
        return true;
    });
}

bool build_vessels(Graph& g, CMap2& m2)
{
    bool okay;
    GraphAttributes gAttribs;
    okay = get_graph_attributes(g, gAttribs);
    if(!okay)
    {
        std::cout << "error build_vessels: get_graph_attributes" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): got graph attributes" << std::endl;

    GraphData gData;
    okay = get_graph_data(g, gData);
    std::cout << gData.intersections.size() << " inters" << std::endl;
    std::cout << gData.branches.size() << " branches" << std::endl;
    for(auto b : gData.branches){
        std::cout << index_of(g, GVertex(b.first)) << " - " << index_of(g, GVertex(b.second)) << std::endl;
    }
    if(!okay)
    {
        std::cout << "error build_vessels: get_graph_data" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): got graph attributes" << std::endl;

    okay = build_connection_interfaces(g, gAttribs, m2);
    if(!okay)
    {
        std::cout << "error build_vessels: add_cmap2_attributes" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): intersection cores built" << std::endl;

    M2Attributes m2Attribs;
    okay = add_cmap2_attributes(m2, m2Attribs);
    if(!okay)
    {
        std::cout << "error build_vessels: add_cmap2_attributes" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): added m2 attributes" << std::endl;

    okay = complete_intersections(g, gAttribs, m2, m2Attribs);
    if(!okay)
    {
        std::cout << "error build_vessels: complete_intersections" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): intersections completed" << std::endl;

    // dump_map(m2);

    okay = create_intersections_frames(g, gAttribs, m2, m2Attribs);
    if(!okay)
    {
        std::cout << "error build_vessels: create_intersections_frames" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): create_intersections_frames completed" << std::endl;
    
    okay = propagate_frames(g, gAttribs, gData, m2);
    if(!okay)
    {
        std::cout << "error build_vessels: propagate_frames" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): propagate_frames completed" << std::endl;
    
    okay = complete_interfaces(g, gAttribs, m2, m2Attribs);
    if(!okay)
    {
        std::cout << "error build_vessels: complete_interfaces" << std::endl;
        return false;
    }
    else
        std::cout << "build_vessels (/): complete_interfaces completed" << std::endl;
      
    return true;    
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_VESSELS_GENERATION_H_
