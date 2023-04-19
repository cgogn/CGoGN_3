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

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/core/types/cmap/graph.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/surface_filtering.h>
#include <cgogn/modeling/ui_modules/skeleton_extractor_module.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
#include <cgogn/rendering/ui_modules/graph_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/core/utils/string.h>

using Graph = cgogn::Graph;
using NonManifold = cgogn::IncidenceGraph;
using Surface = cgogn::CMap2;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string graph_filename, surface_filename;
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " surface_file" << std::endl;
		return 1;
	}
	surface_filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Skeleton extractor");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Graph> mpg(app);
	cgogn::ui::MeshProvider<NonManifold> mpnm(app);
	cgogn::ui::MeshProvider<Surface> mps(app);

	cgogn::ui::SurfaceDifferentialProperties<Surface> sdp(app);
	cgogn::ui::SurfaceModeling<Surface> sm(app);
	cgogn::ui::SurfaceFiltering<Surface> sf(app);

	cgogn::ui::GraphRender<Graph> gr(app);
	cgogn::ui::SurfaceRender<NonManifold> srnm(app);
	cgogn::ui::SurfaceRender<Surface> srs(app);

	cgogn::ui::SkeletonExtractor<Surface, NonManifold, Graph> se(app);

	app.init_modules();

	cgogn::ui::View* v = app.current_view();

	v->link_module(&mpg);
	v->link_module(&mpnm);
	v->link_module(&mps);
	v->link_module(&gr);
	v->link_module(&srnm);
	v->link_module(&srs);

	// load surface
	Surface* s = mps.load_surface_from_file(surface_filename);
	std::cout << "surface loaded" << std::endl;
	if (!s)
	{
		std::cout << "Surface file could not be loaded" << std::endl;
		return 1;
	}

	Surface* clone = mps.clone_mesh(*s);

	srs.set_render_vertices(*v, *s, false);
	srs.set_render_vertices(*v, *clone, false);
	srs.set_render_faces(*v, *clone, false);

	auto surface_vertex_position = cgogn::get_attribute<Vec3, cgogn::mesh_traits<Surface>::Vertex>(*s, "position");
	auto surface_vertex_normal = cgogn::get_or_add_attribute<Vec3, cgogn::mesh_traits<Surface>::Vertex>(*s, "normal");

	sdp.compute_normal(*s, surface_vertex_position.get(), surface_vertex_normal.get());

	return app.launch();
}
