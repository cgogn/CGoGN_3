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

#ifndef CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_
#define CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

template <typename MESH, typename CELL>
class CellMarker;

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;
using Mat3 = geometry::Mat3;

struct GData
{
	std::vector<std::pair<Graph::HalfEdge, Graph::HalfEdge>> branches;
	std::vector<Graph::Vertex> intersections;
};

struct GAttributes
{
	std::shared_ptr<Graph::Attribute<Vec3>> vertex_position;
	std::shared_ptr<Graph::Attribute<Scalar>> vertex_radius;
	std::shared_ptr<Graph::Attribute<Dart>> vertex_contact_surface;
	std::shared_ptr<Graph::Attribute<Dart>> halfedge_contact_surface_face;
	std::shared_ptr<Graph::Attribute<Mat3>> halfedge_frame;
};

struct M2Attributes
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position;
	std::shared_ptr<CMap2::Attribute<Vec3>> volume_center;
	std::shared_ptr<CMap2::Attribute<Vec3>> edge_mid;
	std::shared_ptr<CMap2::Attribute<Dart>> halfedge_volume_connection;
};

bool graph_to_hex(Graph& g, CMap2& m2, CMap3& m3);

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

void index_volume_cells(CMap2& m, CMap2::Volume vol);
Graph::HalfEdge branch_extremity(const Graph& g, Graph::HalfEdge h, CellMarker<Graph, Graph::Edge>& cm);
Dart add_branch_section(CMap3& m3);
void project_on_sphere(Vec3& P, const Vec3& C, Scalar R);
void shift_frame(Mat3& frame, uint32 nb_shifts);

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool get_graph_data(const Graph& g, GData& gData);
bool add_graph_attributes(Graph& g, GAttributes& gAttribs);
bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs);

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::Vertex v);
void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::Vertex v);
void build_contact_surface_3(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v);
void build_contact_surface_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v);

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
bool create_intersection_frame_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v);

bool propagate_frames(const Graph& g, GAttributes& gAttribs, const GData& gData, CMap2& m2);
void propagate_frame_n_1(const Graph& g, GAttributes& gAttribs, Graph::HalfEdge h_from_start);
bool propagate_frame_n_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::HalfEdge h_from_start);

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const Graph& g, const GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/

bool build_branch_sections(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool sew_branch_sections(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);

/*****************************************************************************/
/* utils				                                                     */
/*****************************************************************************/
bool dijkstra_topo(CMap2& m2, std::shared_ptr<CMap2::Attribute<Dart>> previous, std::shared_ptr<CMap2::Attribute<uint>> dist);
bool convex_hull(const Graph& g, GAttributes& gAttribs, Graph::Vertex v, CMap2& hull, M2Attributes& hullAttribs)

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_
