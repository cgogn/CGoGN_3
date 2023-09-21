/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
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

#ifndef CGOGN_GEOMETRY_ALGOS_INVERTION_H_
#define CGOGN_GEOMETRY_ALGOS_INVERTION_H_

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/functions/normal.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
bool edge_has_vertex(MESH& m, typename mesh_traits<MESH>::Edge e, typename mesh_traits<MESH>::Vertex v)
{
	auto iv = incident_vertices(m, e);
	return index_of(m, iv[0]) == index_of(m, v) || index_of(m, iv[1]) == index_of(m, v);
}

template <typename MESH>
bool face_has_vertex(MESH& m, typename mesh_traits<MESH>::Face f, typename mesh_traits<MESH>::Vertex v)
{
	auto iv = incident_vertices(m, f);
	for (auto ifv : iv)
		if (index_of(m, ifv) == index_of(m, v))
			return true;
	return false;
}

template <typename MESH>
bool contractible(MESH& m, typename mesh_traits<MESH>::Edge e,
				  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
				  Vec3 new_vertex_posisiton)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	auto iv = incident_vertices(m, e);
	auto v1 = iv[0];
	auto v2 = iv[1];
	auto v1_index = index_of(m, v1);
	auto v2_index = index_of(m, v2);
	auto ifv1 = incident_faces(m, v1);
	auto ifv2 = incident_faces(m, v2);

	for (Face iface : ifv1)
	{
		if (face_has_vertex(m, iface, v1) && face_has_vertex(m, iface, v2))
			continue;

		std::vector<Vec3> old_face;
		std::vector<Vec3> new_face;
		foreach_incident_vertex(m, iface, [&](Vertex ifacev) {
			/*auto ifaceiv = incident_vertices(m, iface);
			std::set<Vertex> sv(ifaceiv.begin(), ifaceiv.end());
			if (sv.size() != 3)
			{
				std::cout << "face: " << index_of(m, iface) << " is not a triangle" << std::endl;
				for (Vertex ifaceiiv : ifaceiv)
				{
					std::cout << " vertex: " << index_of(m, ifaceiiv) << ", ";
				}
				std::cout << std::endl;
			}*/
			if (index_of(m, ifacev) != v1_index)
			{
				old_face.push_back(value<Vec3>(m, vertex_position, ifacev));
				new_face.push_back(value<Vec3>(m, vertex_position, ifacev));
			}
			return true;
		});
		old_face.push_back(value<Vec3>(m, vertex_position, v1));
		new_face.push_back(new_vertex_posisiton);
		Vec3 old_face_normal = geometry::normal(old_face[0], old_face[1], old_face[2]);
		Vec3 new_face_normal = geometry::normal(new_face[0], new_face[1], new_face[2]);
		if (old_face_normal.dot(new_face_normal) < 0)
			return false;
	}
	for (Face iface : ifv2)
	{
		if (face_has_vertex(m, iface, v1) && face_has_vertex(m, iface, v2))
			continue;
		std::vector<Vec3> old_face;
		std::vector<Vec3> new_face;
		foreach_incident_vertex(m, iface, [&](Vertex ifacev) {
			/*auto ifaceiv= incident_vertices(m, iface);
			std::set<Vertex> sv(ifaceiv.begin(), ifaceiv.end());
			if (sv.size() != 3)
			{
				std::cout << "face: "<<index_of(m, iface)<<" is not a triangle"<<std::endl;
				for (Vertex ifaceiiv : ifaceiv)
				{
					std::cout << " vertex: " << index_of(m, ifaceiiv) << ", ";
				}
				std::cout<< std::endl;
			}*/
			if (index_of(m, ifacev) != v2_index)
			{
				old_face.push_back(value<Vec3>(m, vertex_position, ifacev));
				new_face.push_back(value<Vec3>(m, vertex_position, ifacev));
			}
			return true;
		});
		old_face.push_back(value<Vec3>(m, vertex_position, v2));
		new_face.push_back(new_vertex_posisiton);
		Vec3 old_face_normal = geometry::normal(old_face[0], old_face[1], old_face[2]);
		Vec3 new_face_normal = geometry::normal(new_face[0], new_face[1], new_face[2]);
		if (old_face_normal.dot(new_face_normal) < 0)
			return false;
		// std::cout << "------------------------------------------------------------" << std::endl;
	}
	return true;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_INVERTION_H_
