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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/attributes.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/graph_render/graph_render.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>


#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

using Mesh = cgogn::IncidenceGraph;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " filename" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple graph viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	// cgogn::ui::GraphRender<Mesh> gr(app);
	Mesh* ig0 = mp.load_surface_from_file(filename);

	Mesh ig = Mesh();
	Mesh::Vertex vx0 = cgogn::add_vertex(ig);
	Mesh::Vertex vx1 = cgogn::add_vertex(ig);
	Mesh::Vertex vx2 = cgogn::add_vertex(ig);
	Mesh::Vertex vx3 = cgogn::add_vertex(ig);
	Mesh::Vertex vx4 = cgogn::add_vertex(ig);
	Mesh::Edge e0 = cgogn::add_edge(ig, vx0, vx1);
	Mesh::Edge e1 = cgogn::add_edge(ig, vx2, vx3);
	Mesh::Edge e2 = cgogn::add_edge(ig, vx0, vx3);
	Mesh::Edge e3 = cgogn::add_edge(ig, vx2, vx1);
	Mesh::Edge e4 = cgogn::add_edge(ig, vx0, vx4);
	Mesh::Face f0 = cgogn::add_face(ig, {e0, e1, e2, e3});
	Mesh::Face f1 = cgogn::add_face(ig, {e3, e0, e2, e1});
	Mesh::Face f2 = cgogn::add_face(ig, {e0, e1, e4, e3});

	cgogn::remove_face(ig, f1);
	cgogn::cut_edge(ig, e0);
	foreach_cell(ig, [&](Mesh::Vertex v) -> bool {
		std::cout << v.index_ << " " << (*(ig.vertex_incident_edges_))[v.index_].size() << std::endl;
		return true;
	});

	foreach_cell(ig, [&](Mesh::Edge e) -> bool {
		std::cout << e.index_ << " - " << (*(ig.edge_incident_vertices_))[e.index_].first.index_ << " " << (*(ig.edge_incident_vertices_))[e.index_].second.index_
		<< " - " << (*(ig.edge_incident_faces_))[e.index_].size() << std::endl;
		return true;
	});

	foreach_cell(ig, [&](Mesh::Face f) -> bool {
		std::cout << f.index_ << " " << (*(ig.face_incident_edges_))[f.index_].size() << std::endl;
		for(auto e : (*(ig.face_incident_edges_))[f.index_])
		{
			std::cout << e.index_ << ", ";
		}
		std::cout << std::endl;
		return true;
	});

	// foreach_incident_vertex(ig, f0, [&](Mesh::Vertex v) -> bool {
	// 	std::cout << f0.index_ << " -- " << v.index_ << std::endl;
	// 	return true;
	// });

	// foreach_incident_vertex(ig, e0, [&](Mesh::Vertex v) -> bool {
	// 	std::cout << e0.index_ << " -- " << v.index_ << std::endl;
	// 	return true;
	// });

	// foreach_incident_edge(ig, f0, [&](Mesh::Edge e) -> bool {
	// 	std::cout << f0.index_ << " -+- " << e.index_ << std::endl;
	// 	return true;
	// });

	// std::vector<Mesh::Face> faces = incident_faces(ig, e0);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	// v1->link_module(&mp);
	// v1->link_module(&gr);

	// Mesh* m = mp.load_graph_from_file(filename);
	// if (!m)
	// {
	// 	std::cout << "File could not be loaded" << std::endl;
	// 	return 1;
	// }

	// std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	// std::shared_ptr<Attribute<Scalar>> vertex_radius = cgogn::get_attribute<Scalar, Vertex>(*m, "radius");
	// mp.set_mesh_bb_vertex_position(m, vertex_position);
	// gr.set_vertex_position(*v1, *m, vertex_position);
	// gr.set_vertex_radius(*v1, *m, vertex_radius);

	return app.launch();
}
