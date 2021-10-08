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

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/modeling/algos/convex_hull.h>

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;

int main(int argc, char** argv)
{
	using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
	using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
	using Face = typename cgogn::mesh_traits<Mesh>::Face;

	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Sphere partition");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);

	std::string filename;
	if (argc > 1)
		filename = std::string(argv[1]);
	Mesh* m = mp.load_surface_from_file(filename);

	Mesh* hull = mp.add_mesh("hull");
	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::add_attribute<Vec3, Vertex>(*hull, "position");

	if (m)
	{
		std::shared_ptr<Attribute<Vec3>> m_vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
		std::vector<Vec3> points;
		points.reserve(mp.mesh_data(*m).nb_cells<Vertex>());
		for (const Vec3& p : *m_vertex_position)
			points.push_back(p);
		cgogn::modeling::convex_hull(points, *hull, vertex_position.get());

		mp.set_mesh_bb_vertex_position(*m, m_vertex_position);
		sr.set_vertex_position(*v1, *m, m_vertex_position);
	}
	else
		cgogn::modeling::convex_hull({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, *hull, vertex_position.get());

	mp.emit_connectivity_changed(*hull);
	mp.emit_attribute_changed(*hull, vertex_position.get());

	mp.set_mesh_bb_vertex_position(*hull, vertex_position);
	sr.set_vertex_position(*v1, *hull, vertex_position);

	return app.launch();
}
