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

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/Gregory_render/gregory_render_pc.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
using Mesh = cgogn::CMap2;

//#include <cgogn/ui/modules/surface_selection/surface_selection.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;

int main(int argc, char** argv)
{
	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	std::string filename;
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/simple/cube.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Gregory");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::GregoryRenderPC<Mesh> gr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&gr);

	// cgogn::ui::View* v2 = app.add_view();
	// v2->link_module(&mp);
	// v2->link_module(&gr);

	Mesh* m = mp.load_surface_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

	mp.set_mesh_bb_vertex_position(m, vertex_position);
	sr.set_vertex_position(*v1, *m, vertex_position);

	return app.launch();
}
