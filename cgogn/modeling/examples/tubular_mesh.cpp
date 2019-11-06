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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/cmap/cmap_info.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/graph_render/graph_render.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
// #include <cgogn/ui/modules/volume_render/volume_render.h>

#include <cgogn/modeling/algos/vessels_generation.h>


#include <cgogn/modeling/algos/vessels_generation.h>


using namespace cgogn::numerics;

using Graph = cgogn::Graph;
using Surface = cgogn::CMap2;
using Volume = cgogn::CMap3;

template <typename T>
using GraphAttribute = typename cgogn::mesh_traits<Graph>::Attribute<T>;
template <typename T>
using SurfaceAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;
template <typename T>
using VolumeAttribute = typename cgogn::mesh_traits<Volume>::Attribute<T>;

using GraphVertex = typename cgogn::mesh_traits<Graph>::Vertex;
using SurfaceVertex = typename cgogn::mesh_traits<Surface>::Vertex;
using VolumeVertex = typename cgogn::mesh_traits<Volume>::Vertex;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

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
	cgogn::ui::SurfaceRender<Volume> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mpg);
	v1->link_module(&mps);
	v1->link_module(&mpv);

	v1->link_module(&gr);
	v1->link_module(&sr);
	v1->link_module(&vr);

	Graph* g = mpg.load_graph_from_file(filename);
	if (!g)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	Surface* s = mps.add_mesh("contact");
	Volume* v = mpv.add_mesh("tubes");
	Graph* g2 = mpg.add_mesh("vertices");
	std::shared_ptr<Graph::Attribute<cgogn::geometry::Vec3>> vertex_position = cgogn::get_attribute<cgogn::geometry::Vec3, GraphVertex>(*g, "position");
	// std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");

	// mp.set_mesh_bb_vertex_position(m, vertex_position);

	// sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());
	
	// gr.set_vertex_position(*g, vertex_position);
	// sr.set_vertex_normal(*m, vertex_normal);
	// auto vertices = 
	// auto edges = cgogn::incident_edges()

	// cgogn::index_cells<Graph::Edge>(*g);
	// cgogn::foreach_cell(*g, [&](Graph::Vertex v) ->bool 
	// {
	// 	std::cout << cgogn::index_of(*g, v) << "-" ;
	// 	auto edges = cgogn::incident_edges(*g, v);
	// 	for(auto e : edges)
	// 		std::cout << cgogn::index_of(*g, e) << " ";
	// 	std::cout << std::endl;
	// 	return true;
	// });


	if(cgogn::modeling::build_vessels(*g, *s, *v, *g2))
	{
		mpg.emit_connectivity_changed(g);
		mpg.emit_attribute_changed(g, vertex_position.get());

	std::shared_ptr<Graph::Attribute<cgogn::geometry::Vec3>> vertex_position2 = cgogn::get_attribute<cgogn::geometry::Vec3, GraphVertex>(*g2, "position");
		mpg.set_mesh_bb_vertex_position(g2, vertex_position2);

		mpg.emit_connectivity_changed(g2);
		mpg.emit_attribute_changed(g2, vertex_position2.get());

		std::shared_ptr<Surface::Attribute<cgogn::geometry::Vec3>> vertex_position_s = cgogn::get_attribute<cgogn::geometry::Vec3, SurfaceVertex>(*s, "position");
	// mp.set_mesh_bb_vertex_position(m, vertex_position);
		mps.emit_connectivity_changed(s);
		mps.emit_attribute_changed(s, vertex_position_s.get());

		std::shared_ptr<Volume::Attribute<cgogn::geometry::Vec3>> vertex_position_v = cgogn::get_attribute<cgogn::geometry::Vec3, VolumeVertex>(*v, "position");
		mpv.emit_connectivity_changed(v);
		mpv.emit_attribute_changed(v, vertex_position_v.get());
	}
	// dump_map(*s);

	return app.launch();
}
