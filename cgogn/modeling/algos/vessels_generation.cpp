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

#include <cgogn/modeling/algos/vessels_generation.h>

namespace cgogn
{

namespace modeling
{


using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;
// using GVertex = Graph::Vertex;
// using GEdge = Graph::Edge;
// using M2Vertex = CMap2::Vertex;
// using M2Face = CMap2::Face;

bool get_graph_attributes(const Graph& g, GraphAttributes& gAttribs)
{
    gAttribs.V_pos = cgogn::get_attribute<Vec3, Graph::Vertex>(g, "position");
    if(!gAttribs.V_pos)
    {
        std::cout << "The graph has no position atttribute" << std::endl;
        return false;
    }

    gAttribs.V_rad = cgogn::get_attribute<Scalar, Graph::Vertex>(g, "radius");
    if(!gAttribs.V_rad)
    {
        std::cout << "The graph has no radius atttribute" << std::endl;
        return false;
    }

    return true;
}

bool clean_graph(Graph& g)
{
    return true;
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
        std::cout << "build_vessels: got graph attributes" << std::endl;

    // std::shared_ptr<CMap2::Attribute<Vec3>> M2_v_pos = cgogn::add_attribute<Vec3, CMap2::Vertex>(m2, "position");
    // std::shared_ptr<Graph::Attribute<Vec3>> G_v_pos = cgogn::get_attribute<Vec3, Graph::Vertex>(g, "position");
    // CMap2::Face f = add_face(m2, 4, true);
    
    // Vec3 pos = {0.1, 0.1, 0.1};
    // std::vector<CMap2::Vertex> verts = incident_vertices(m2, f);
    // uint32 i = 0;
    // foreach_cell(g, [&](Graph::Vertex v) -> bool
    // {
    //     value<Vec3>(m2, M2_v_pos, verts[i++]) = value<Vec3>(g, G_v_pos, v);
    //     return true;
    // });

	// CellCache<Graph> cache(g);
	// cache.build<Graph::Edge>();
    // foreach_cell(cache, [&](Graph::Edge ge) -> bool
    // {
    //     // value<Vec3>(m2, M2_v_pos, verts[i++]) = value<Vec3>(g, G_v_pos, v);
    //     auto verts = incident_vertices(g, ge);
    //     // std::cout <<
    //     Vec3 midPoint = 0.5 * (value<Vec3>(g, G_v_pos, verts[0]) + value<Vec3>(g, G_v_pos, verts[1])); 
    //     Graph::Vertex v = cut_edge(g, ge);
    //     std::cout << index_of(g, v) << std::endl
    //     << midPoint[0] << "/" << midPoint[1] << "/" << midPoint[2] << std::endl;
    //     value<Vec3>(g, G_v_pos, v) = midPoint;
    //     return true;
    // });

    // foreach_cell(g, [&](Graph::Vertex gv) -> bool
    // {
    //     Vec3 pos = value<Vec3>(g, G_v_pos, gv);
    //     std::cout << index_of(g, gv) << " -- "
    //     << pos[0] << "/" << pos[1] << "/" << pos[2] << std::endl;
    //     return true;
    // });


    // for(uint32 i = 0; i < verts.size(); ++i)
    // {
    //     value<Vec3>(m2, vertex_pos, verts[i]) = i * pos; 
    // }
    // foreach_incident_vertex(m2, f, [&](CMap2::Vertex v) -> bool
    // {
    //     std::cout << m2.embedding(v) << std::endl;

    //     return true;
    // });

    return true;    
}

} // namespace modeling

} // namespace cgogn