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

#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>
#include <cgogn/modeling/algos/incidenceGraph_to_hex.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

using IGraph = cgogn::IncidenceGraph;
using Surface = cgogn::CMap2;
using Volume = cgogn::CMap3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<IGraph>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<IGraph>::Vertex;

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
	app.set_window_title("Simple incidence graph viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<IGraph> mpig(app);
	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Volume> mpv(app);

	cgogn::ui::SurfaceRender<IGraph> gr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::SurfaceRender<Volume> vr(app);

	IGraph* ig = mpig.load_surface_from_file(filename);
	std::cout << filename << " : loaded" << std::endl;

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*ig, "position");

	app.init_modules();
	std::cout << "init_modules : done" << std::endl;

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mpig);
	v1->link_module(&mps);
	v1->link_module(&mpv);
	v1->link_module(&gr);
	v1->link_module(&sr);
	v1->link_module(&vr);

	std::cout << "link_modules : done" << std::endl;

	mpig.set_mesh_bb_vertex_position(*ig, vertex_position);
	gr.set_vertex_position(*v1, *ig, vertex_position);

	std::cout << "display : done" << std::endl;


	Surface* m2 = mps.add_mesh("surface");
	Volume* m3 = mpv.add_mesh("volumes");
	cgogn::modeling::incidenceGraph_to_hex(*ig, *m2, *m3);

	std::shared_ptr<Attribute<Vec3>> surface_vertex_position = cgogn::get_attribute<Vec3, Surface::Vertex>(*m2, "position");

	mps.set_mesh_bb_vertex_position(*m2, surface_vertex_position);
	mps.emit_connectivity_changed(*m2);

	sr.set_vertex_position(*v1, *m2, surface_vertex_position);
	sr.set_render_vertices(*v1, *m2, false);
	sr.set_render_faces(*v1, *m2, false);

	mpv.set_mesh_bb_vertex_position(*m3, surface_vertex_position);
	mpv.emit_connectivity_changed(*m3);

	return app.launch();
}
