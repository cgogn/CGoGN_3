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

#include <cgogn/core/utils/string.h>

#include <cgogn/geometry/algos/distance.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/algos/filtering.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/picking.h>

#include <cgogn/modeling/algos/graph_resampling.h>
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
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string graph_filename, surface_filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " graph_file [enclosing_surface_file]" << std::endl;
		return 1;
	}
	if (argc >= 2)
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

	cgogn::ui::GraphRender<Graph> gr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::VolumeRender<Volume> vr(app);

	app.init_modules();

	cgogn::ui::View* v = app.current_view();

	v->link_module(&mpg);
	v->link_module(&mps);
	v->link_module(&mpv);

	v->link_module(&gr);
	v->link_module(&sr);
	v->link_module(&vr);

	// load graph
	Graph* g = mpg.load_graph_from_file(graph_filename);
	if (!g)
	{
		std::cout << "Graph file could not be loaded" << std::endl;
		return 1;
	}
	auto g_vertex_position = cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Graph>::Vertex>(*g, "position");
	auto g_vertex_radius = cgogn::get_attribute<Scalar, typename cgogn::mesh_traits<Graph>::Vertex>(*g, "radius");
	if (!g_vertex_radius)
	{
		g_vertex_radius = cgogn::add_attribute<Scalar, typename cgogn::mesh_traits<Graph>::Vertex>(*g, "radius");
		Scalar l = cgogn::geometry::mean_edge_length(*g, g_vertex_position.get());
		g_vertex_radius->fill(l);
	}

	// create resampled graph
	Graph* rg = mpg.add_mesh(cgogn::filename_from_path(graph_filename) + "_resampled");
	auto rg_vertex_position = cgogn::add_attribute<Vec3, typename cgogn::mesh_traits<Graph>::Vertex>(*rg, "position");
	auto rg_vertex_radius = cgogn::add_attribute<Scalar, typename cgogn::mesh_traits<Graph>::Vertex>(*rg, "radius");

	// load enclosing surface
	Surface* enclosing = mps.load_surface_from_file(surface_filename);
	if (enclosing)
	{
		auto enclosing_vertex_position =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*enclosing, "position");
		cgogn::modeling::compute_graph_radius_from_surface(*g, g_vertex_position.get(), g_vertex_radius.get(),
														   *enclosing, enclosing_vertex_position.get());
	}

	// resample graph
	cgogn::modeling::resample_graph(*g, g_vertex_position.get(), g_vertex_radius.get(), *rg, rg_vertex_position.get(),
									rg_vertex_radius.get());
	Scalar min_radius = std::numeric_limits<Scalar>::max();
	for (Scalar r : *rg_vertex_radius)
		if (r < min_radius)
			min_radius = r;
	rg_vertex_radius->fill(min_radius);

	mpg.emit_connectivity_changed(g);
	mpg.emit_attribute_changed(g, g_vertex_position.get());
	mpg.emit_attribute_changed(g, g_vertex_radius.get());
	mpg.emit_connectivity_changed(rg);

	// create contact surface
	Surface* contact = mps.add_mesh("contact");
	// create hex mesh
	Volume* hex = mpv.add_mesh("hex");

	bool hex_built = cgogn::modeling::graph_to_hex(*rg, *contact, *hex);
	mps.emit_connectivity_changed(contact);
	mpv.emit_connectivity_changed(hex);

	if (hex_built && enclosing)
	{
		Surface* skin = mps.add_mesh("skin");
		auto skin_vertex_position =
			cgogn::add_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*skin, "position");
		auto skin_filtered_vertex_position =
			cgogn::add_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*skin, "__filtered_position");
		auto skin_vertex_normal =
			cgogn::add_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*skin, "normal");
		auto skin_vertex_laplacian =
			cgogn::add_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*skin, "laplacian");
		auto skin_vertex_hex_vertex =
			cgogn::add_attribute<typename cgogn::mesh_traits<Volume>::Vertex,
								 typename cgogn::mesh_traits<Surface>::Vertex>(*skin, "hex_vertex");

		auto hex_vertex_position =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Volume>::Vertex>(*hex, "position");

		cgogn::modeling::extract_volume_surface(*hex, hex_vertex_position.get(), *skin, skin_vertex_position.get(),
												skin_vertex_hex_vertex.get());

		cgogn::geometry::apply_ear_triangulation(*skin, skin_vertex_position.get());

		using SelectedFace = std::tuple<typename cgogn::mesh_traits<Surface>::Face, Vec3, Scalar>;

		auto enclosing_vertex_position =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*enclosing, "position");

		cgogn::geometry::compute_normal(*skin, skin_vertex_position.get(), skin_vertex_normal.get());

		cgogn::parallel_foreach_cell(*skin, [&](typename cgogn::mesh_traits<Surface>::Vertex v) -> bool {
			Vec3& p = cgogn::value<Vec3>(*skin, skin_vertex_position, v);
			const Vec3& n = cgogn::value<Vec3>(*skin, skin_vertex_normal, v);
			std::vector<SelectedFace> selectedfaces =
				cgogn::geometry::internal::picking(*enclosing, enclosing_vertex_position.get(), p, p + n);
			Vec3 pos = selectedfaces.size() > 0 ? std::get<1>(selectedfaces[0]) : p;
			cgogn::value<Vec3>(
				*hex, hex_vertex_position,
				cgogn::value<typename cgogn::mesh_traits<Volume>::Vertex>(*skin, skin_vertex_hex_vertex, v)) = pos;
			p = pos;
			return true;
		});

		for (cgogn::uint32 i = 0; i < 3; ++i)
		{
			// cgogn::geometry::filter_average<Vec3>(*skin, skin_vertex_position.get(),
			// 									  skin_filtered_vertex_position.get());
			// skin_vertex_position->swap(skin_filtered_vertex_position.get());
			cgogn::geometry::compute_laplacian(*skin, skin_vertex_position.get(), skin_vertex_laplacian.get());
			parallel_foreach_cell(*skin, [&](typename cgogn::mesh_traits<Surface>::Vertex v) -> bool {
				cgogn::value<Vec3>(*skin, skin_vertex_position, v) +=
					0.1 * cgogn::value<Vec3>(*skin, skin_vertex_laplacian, v);
				return true;
			});
			cgogn::parallel_foreach_cell(*skin, [&](typename cgogn::mesh_traits<Surface>::Vertex v) -> bool {
				Vec3 p = cgogn::geometry::closest_point_on_surface(*enclosing, enclosing_vertex_position.get(),
																   cgogn::value<Vec3>(*skin, skin_vertex_position, v));
				cgogn::value<Vec3>(*skin, skin_vertex_position, v) = p;
				cgogn::value<Vec3>(
					*hex, hex_vertex_position,
					cgogn::value<typename cgogn::mesh_traits<Volume>::Vertex>(*skin, skin_vertex_hex_vertex, v)) = p;
				return true;
			});
		}

		mps.emit_connectivity_changed(skin);
	}

	// set bounding box & render attributes

	mpg.set_mesh_bb_vertex_position(g, g_vertex_position);
	mpg.set_mesh_bb_vertex_position(rg, rg_vertex_position);
	gr.set_vertex_position(*v, *g, g_vertex_position);
	gr.set_vertex_position(*v, *rg, rg_vertex_position);

	if (enclosing)
	{
		auto enclosing_vertex_position =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Surface>::Vertex>(*enclosing, "position");
		mps.set_mesh_bb_vertex_position(enclosing, enclosing_vertex_position);
		sr.set_vertex_position(*v, *enclosing, enclosing_vertex_position);
	}

	if (hex_built)
	{
		auto hex_vertex_position =
			cgogn::get_attribute<Vec3, typename cgogn::mesh_traits<Volume>::Vertex>(*hex, "position");
		mpv.set_mesh_bb_vertex_position(hex, hex_vertex_position);
		vr.set_vertex_position(*v, *hex, hex_vertex_position);
	}

	return app.launch();
}
