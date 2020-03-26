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

#include <cgogn/ui/modules/graph_render/graph_render.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>

#include <cgogn/modeling/algos/graph_to_hex.h>

using Graph = cgogn::Graph;
using Surface = cgogn::CMap2;
using Volume = cgogn::CMap3;

template <typename T>
using GraphAttribute = typename cgogn::mesh_traits<Graph>::Attribute<T>;

template <typename T>
using SurfaceAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;

template <typename T>
using VolumeAttribute = typename cgogn::mesh_traits<Volume>::Attribute<T>;

using Vec3 = cgogn::geometry::Vec3;

int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " graph_file" << std::endl;
		return 1;
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Tubular mesh");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Graph> mpg(app);
	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Volume> mpv(app);

	cgogn::ui::GraphRender<Graph> gr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	//cgogn::ui::VolumeRender<Volume> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mpg);
	v1->link_module(&mps);
	v1->link_module(&mpv);

	v1->link_module(&gr);
	v1->link_module(&sr);
	//v1->link_module(&vr);

	Graph* g = mpg.load_graph_from_file(filename);
	if (!g)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}
	Surface* s = mps.add_mesh("contact");
	Volume* v = mpv.add_mesh("tubes");

	if (cgogn::modeling::graph_to_hex(*g, *s, *v))
	{
		std::shared_ptr<GraphAttribute<Vec3>> vertex_position_g =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Graph>::Vertex>(*g, "position");
		mpg.emit_connectivity_changed(g);
		mpg.emit_attribute_changed(g, vertex_position_g.get());

		std::shared_ptr<SurfaceAttribute<Vec3>> vertex_position_s =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*s, "position");
		mps.emit_connectivity_changed(s);
		mps.emit_attribute_changed(s, vertex_position_s.get());

		//std::shared_ptr<VolumeAttribute<Vec3>> vertex_position_v =
		//	cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Volume>::Vertex>(*v, "position");
		//vr.set_vertex_position(*v1, *v, vertex_position_v);

		//mpv.emit_connectivity_changed(v);
		//mpv.emit_attribute_changed(v, vertex_position_v.get());

	}

	return app.launch();
}
