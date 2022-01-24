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
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/vector_per_vertex_render.h>
// #include <cgogn/geometry/ui_modules/surface_selection.h>

#include <cgogn/core/types/mesh_views/cell_filter.h>
#include <cgogn/modeling/algos/subdivision.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

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

	std::string filename;
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Modeling");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::VectorPerVertexRender<Mesh> vpvr(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);
	cgogn::ui::SurfaceModeling<Mesh> sm(app);
	// cgogn::ui::SurfaceSelection<Mesh> ss(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&vpvr);
	// v1->link_module(&ss);

	Mesh* m = mp.load_surface_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");

	mp.set_mesh_bb_vertex_position(*m, vertex_position);

	sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());

	// cgogn::CellFilter<Mesh> cf(*m);
	// cf.set_filter<Edge>([](Edge) -> bool { return true; });
	// cf.set_filter<Face>([&](Face f) -> bool {
	// 	std::vector<Vertex> iv = incident_vertices(*m, f);
	// 	return cgogn::value<Vec3>(*m, vertex_position, iv[0])[0] < 0;
	// });

	// cgogn::CellCache<Mesh> cc(*m);
	// cc.build<Edge>();
	// cc.build<Face>([&](Face f) -> bool {
	// 	std::vector<Vertex> iv = incident_vertices(*m, f);
	// 	return cgogn::value<Vec3>(*m, vertex_position, iv[0])[0] < 0;
	// });

	// cgogn::modeling::quadrangulate_all_faces(
	// 	cc,
	// 	[&](Vertex v) {
	// 		std::vector<Vertex> av = adjacent_vertices_through_edge(*m, v);
	// 		cgogn::value<Vec3>(*m, vertex_position, v) =
	// 			0.5 * (cgogn::value<Vec3>(*m, vertex_position, av[0]) + cgogn::value<Vec3>(*m, vertex_position, av[1]));
	// 	},
	// 	[&](Vertex v) {
	// 		Vec3 center;
	// 		center.setZero();
	// 		uint32 count = 0;
	// 		foreach_adjacent_vertex_through_edge(*m, v, [&](Vertex av) -> bool {
	// 			center += cgogn::value<Vec3>(*m, vertex_position, av);
	// 			++count;
	// 			return true;
	// 		});
	// 		center /= Scalar(count);
	// 		cgogn::value<Vec3>(*m, vertex_position, v) = center;
	// 	});

	// mp.emit_connectivity_changed(m);
	// mp.emit_attribute_changed(m, vertex_position.get());

	sr.set_vertex_position(*v1, *m, vertex_position);
	sr.set_vertex_normal(*v1, *m, vertex_normal);

	return app.launch();
}
