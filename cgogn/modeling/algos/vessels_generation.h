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
// using GVertex = Graph::Vertex;
// using GEdge = Graph::Edge;
// using M2Vertex = CMap2::Vertex;
// using M2Face = CMap2::Face;

struct GraphAttributes
{
    // typename mesh_traits<Graph>::template Graph::Attribute<Vec3>* vertex_position;

};

bool get_graph_attributes(const Graph& g, GraphAttributes& g_attribs)
{

}

bool clean_graph(Graph& g)
{
    return true;
}



void build_vessels(Graph& g, CMap2& m2)
{
    std::shared_ptr<CMap2::Attribute<Vec3>> M2_v_pos = cgogn::add_attribute<Vec3, CMap2::Vertex>(m2, "position");
    std::shared_ptr<Graph::Attribute<Vec3>> G_v_pos = cgogn::get_attribute<Vec3, Graph::Vertex>(g, "position");
    CMap2::Face f = add_face(m2, 3, true);
    
    Vec3 pos = {0.1, 0.1, 0.1};
    std::vector<CMap2::Vertex> verts = incident_vertices(m2, f);
    uint32 i = 0;
    foreach_cell(g, [&](Graph::Vertex v) -> bool
    {
        value<Vec3>(m2, M2_v_pos, verts[i++]) = value<Vec3>(g, G_v_pos, v);
        return true;
    });

    // for(uint32 i = 0; i < verts.size(); ++i)
    // {
    //     value<Vec3>(m2, vertex_pos, verts[i]) = i * pos; 
    // }
    // foreach_incident_vertex(m2, f, [&](CMap2::Vertex v) -> bool
    // {
    //     std::cout << m2.embedding(v) << std::endl;

    //     return true;
    // });

    return;    
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_VESSELS_GENERATION_H_
