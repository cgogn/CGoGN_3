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
    m2Attribs.V_pos = cgogn::add_attribute<Vec3, CMap2::Vertex>(m2, "position");
    if(!m2Attribs.V_pos)
    {
        std::cout << "Failed to add position attribute to map2" << std::endl;
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
                build_interface_core(g, gAttribs, m2, gv);
                break;
            default:
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
                //complete_intersection_n
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
        Dart d0 = value<Dart>(g, gAttribs.m2_interface, GDart(d));
        Dart d1 = m2.phi<11>(d0);

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



bool complete_interface_1(const Graph& g, const GraphAttributes& gAttribs, 
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

    okay = create_intersections_frames(g, gAttribs, m2, m2Attribs);

    okay = complete_interfaces(g, gAttribs, m2, m2Attribs);

    return true;    
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_VESSELS_GENERATION_H_
