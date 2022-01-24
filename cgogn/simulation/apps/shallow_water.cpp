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
#include <cgogn/geometry/ui_modules/surface_selection.h>
#include <cgogn/modeling/ui_modules/surface_deformation.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/vector_per_vertex_render.h>
#include <cgogn/simulation/ui_modules/shallow_water.h>

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;

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
	app.set_window_title("Shallow water");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::VectorPerVertexRender<Mesh> vpvr(app);
	cgogn::ui::SurfaceSelection<Mesh> ss(app);
	cgogn::ui::SurfaceDeformation<Mesh> sd(app);
	cgogn::ui::ShallowWater<Mesh> sw(app);

	cgogn::ui::View* v1 = app.current_view();
	cgogn::ui::View* v2 = app.add_view();

	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&vpvr);

	v2->link_module(&mp);
	v2->link_module(&sr);
	v2->link_module(&ss);
	v2->link_module(&sd);

	app.init_modules();

	Mesh* domain = mp.load_surface_from_file(filename);
	if (!domain)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	sw.set_domain(domain);

	return app.launch();
}
