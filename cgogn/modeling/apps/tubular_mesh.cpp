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

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap3.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
// #include <cgogn/geometry/ui_modules/registration.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/volume_selection.h>
#include <cgogn/modeling/ui_modules/tubular_mesh_module.h>
#include <cgogn/modeling/ui_modules/volume_deformation.h>
// #include <cgogn/rendering/ui_modules/graph_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/vector_per_vertex_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>

#include <cgogn/core/utils/string.h>

// using Graph = cgogn::Graph;
using Graph = cgogn::IncidenceGraph;
using Surface = cgogn::CMap2;
using Volume = cgogn::CMap3;

template <typename T>
using GraphAttribute = typename cgogn::mesh_traits<Graph>::Attribute<T>;
template <typename T>
using SurfaceAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;
template <typename T>
using VolumeAttribute = typename cgogn::mesh_traits<Volume>::Attribute<T>;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string graph_filename, surface_filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " graph_file [enclosing_surface_file]" << std::endl;
		return 1;
	}
	graph_filename = std::string(argv[1]);
	if (argc >= 3)
		surface_filename = std::string(argv[2]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Tubular mesh");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Graph> mpg(app);
	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Volume> mpv(app);

	cgogn::ui::VolumeSelection<Volume> vs(app);
	cgogn::ui::VolumeDeformation<Volume> vd(app);

	cgogn::ui::SurfaceDifferentialProperties<Surface> sdp(app);

	// cgogn::ui::Registration<Surface> r(app);

	cgogn::ui::SurfaceRender<Graph> gr(app);
	// cgogn::ui::GraphRender<Graph> gr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::VectorPerVertexRender<Surface> svpvr(app);
	cgogn::ui::VectorPerVertexRender<Graph> gvpvr(app);
	cgogn::ui::VolumeRender<Volume> vr(app);

	cgogn::ui::TubularMesh<Graph, Surface, Volume> tm(app);

	app.init_modules();

	cgogn::ui::View* v = app.current_view();

	v->link_module(&mpg);
	v->link_module(&mps);
	v->link_module(&mpv);

	v->link_module(&vs);
	v->link_module(&vd);

	v->link_module(&gr);
	v->link_module(&sr);
	v->link_module(&svpvr);
	v->link_module(&gvpvr);
	v->link_module(&vr);

	v->link_module(&tm);

	// load graph
	std::cout << "load graph" << std::endl;
	Graph* g = mpg.load_graph_from_file(graph_filename);
	std::cout << "graph loaded" << std::endl;
	if (!g)
	{
		std::cout << "Graph file could not be loaded" << std::endl;
		return 1;
	}

	// load surface
	if (argc >= 3)
	{
		Surface* s = mps.load_surface_from_file(surface_filename, false);
		std::cout << "surface loaded" << std::endl;
		if (!s)
		{
			std::cout << "Surface file could not be loaded" << std::endl;
			return 1;
		}

		sr.set_render_vertices(*v, *s, false);
		sr.set_render_faces(*v, *s, false);

		auto graph_vertex_position = cgogn::get_attribute<Vec3, cgogn::mesh_traits<Graph>::Vertex>(*g, "position");
		auto graph_vertex_radius = cgogn::get_attribute<Scalar, cgogn::mesh_traits<Graph>::Vertex>(*g, "radius");
		if (!graph_vertex_radius)
			graph_vertex_radius = cgogn::add_attribute<Scalar, cgogn::mesh_traits<Graph>::Vertex>(*g, "radius");

		tm.set_current_graph(g);
		tm.set_current_graph_vertex_position(graph_vertex_position);
		tm.set_current_graph_vertex_radius(graph_vertex_radius);

		auto surface_vertex_position = cgogn::get_attribute<Vec3, cgogn::mesh_traits<Surface>::Vertex>(*s, "position");

		tm.set_current_surface(s);
		tm.set_current_surface_vertex_position(surface_vertex_position);

		// tm.extend_graph_extremities();
		// tm.extend_graph_extremities();
		// tm.init_graph_radius_from_surface();
		// Graph* resampled_graph = tm.resample_graph();

		// auto resampled_graph_vertex_position =
		// 	cgogn::get_attribute<Vec3, cgogn::mesh_traits<Graph>::Vertex>(*g, "position");
		// auto resampled_graph_vertex_radius =
		// 	cgogn::get_attribute<Scalar, cgogn::mesh_traits<Graph>::Vertex>(*g, "radius");

		// tm.set_current_graph(resampled_graph);
		// tm.set_current_graph_vertex_position(resampled_graph_vertex_position);
		// tm.set_current_graph_vertex_radius(resampled_graph_vertex_radius);

		// *********** //

		// auto start_timer = std::chrono::high_resolution_clock::now();

		// Volume* h = tm.build_hex_mesh();

		// auto end_timer = std::chrono::high_resolution_clock::now();
		// std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
		// std::cout << "hex mesh generation in " << elapsed_seconds.count() << " seconds" << std::endl;

		// start_timer = end_timer;

		// for (int i = 0; i < 30; ++i)
		// 	tm.regularize_surface_vertices(5.0);
		// tm.add_volume_padding(true);
		// for (int i = 0; i < 5; ++i)
		// 	tm.relocate_interior_vertices();
		// for (int i = 0; i < 15; ++i)
		// {
		// 	tm.refresh_edge_target_length_ = true;
		// 	tm.optimize_volume_vertices(1.0, false);
		// }
		// tm.subdivide_volume();
		// for (int i = 0; i < 20; ++i)
		// 	tm.optimize_volume_vertices(1.0, false);
		// tm.optimize_volume_vertices(1.0, true);
		// // for (int i = 0; i < 10; ++i)
		// // 	tm.optimize_volume_vertices(12.0, true);

		// end_timer = std::chrono::high_resolution_clock::now();
		// elapsed_seconds = end_timer - start_timer;
		// std::cout << "regul + subdiv + optim in " << elapsed_seconds.count() << " seconds" << std::endl;

		// auto volume_vertex_position = cgogn::get_attribute<Vec3, cgogn::mesh_traits<Volume>::Vertex>(*h, "position");
		// vr.set_vertex_position(*v, *h, volume_vertex_position);

		// *********** //
	}

	// load volume
	if (argc >= 4)
	{
		std::string volume_filename = std::string(argv[3]);
		Volume* h = mpv.load_volume_from_file(volume_filename);
		if (!h)
		{
			std::cout << "Volume file could not be loaded" << std::endl;
			return 1;
		}

		tm.set_current_volume(h);
		auto volume_vertex_position = cgogn::get_attribute<Vec3, cgogn::mesh_traits<Volume>::Vertex>(*h, "position");
		vr.set_vertex_position(*v, *h, volume_vertex_position);
	}

	return app.launch();
}
