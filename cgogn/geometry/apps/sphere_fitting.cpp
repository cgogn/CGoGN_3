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
#include <cgogn/core/types/maps/cmap/cmap2.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/sphere_fitting.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

// #include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

using Surface = cgogn::CMap2;
using Points = cgogn::CMap0;
using NonManifold = cgogn::IncidenceGraph;

template <typename T>
using SAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;
using SVertex = typename cgogn::mesh_traits<Surface>::Vertex;

template <typename T>
using PAttribute = typename cgogn::mesh_traits<Points>::Attribute<T>;
using PVertex = typename cgogn::mesh_traits<Points>::Vertex;

template <typename T>
using NMAttribute = typename cgogn::mesh_traits<NonManifold>::Attribute<T>;
using NMVertex = typename cgogn::mesh_traits<NonManifold>::Vertex;

using cgogn::geometry::Scalar;
using cgogn::geometry::Vec3;
using cgogn::geometry::Vec4;

int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Sphere fitting");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Points> mpp(app);
	cgogn::ui::MeshProvider<NonManifold> mpnm(app);

	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::SurfaceRender<NonManifold> srnm(app);
	cgogn::ui::PointCloudRender<Points> pcr(app);

	cgogn::ui::SurfaceModeling<Surface> sm(app);
	cgogn::ui::SphereFitting<Surface, Points, NonManifold> sf(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mps);
	v1->link_module(&mpp);
	v1->link_module(&sr);
	v1->link_module(&srnm);
	v1->link_module(&pcr);
	v1->link_module(&sf);

	Surface* s = mps.load_surface_from_file(filename);
	if (!s)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	auto s_vertex_position = cgogn::get_attribute<Vec3, SVertex>(*s, "position");

	mps.set_mesh_bb_vertex_position(*s, s_vertex_position);

	// cgogn::modeling::pliant_remeshing(*s, s_vertex_position, 1.0, true, true);
	// mps.emit_connectivity_changed(*s);
	// mps.emit_attribute_changed(*s, s_vertex_position.get());

	sf.set_selected_surface(*s);
	sf.set_surface_vertex_position(*s, s_vertex_position);
	sf.init_surface_data(*s);
	sf.init_spheres(*s);

	auto s_vertex_normal = cgogn::get_attribute<Vec3, SVertex>(*s, "normal");

	sr.set_vertex_position(*v1, *s, s_vertex_position);
	sr.set_vertex_normal(*v1, *s, s_vertex_normal);
	sr.set_render_vertices(*v1, *s, false);
	sr.set_ghost_mode(*v1, *s, true);
	sr.set_render_edges(*v1, *s, false);

	Points* spheres = mpp.mesh(mps.mesh_name(*s) + "_spheres");
	auto p_vertex_position = cgogn::get_attribute<Vec3, PVertex>(*spheres, "position");
	auto p_vertex_radius = cgogn::get_attribute<Scalar, PVertex>(*spheres, "radius");
	auto p_vertex_color = cgogn::get_attribute<Vec4, PVertex>(*spheres, "color");

	pcr.set_vertex_position(*v1, *spheres, p_vertex_position);
	pcr.set_vertex_radius(*v1, *spheres, p_vertex_radius);
	pcr.set_vertex_color(*v1, *spheres, p_vertex_color);
	pcr.set_vertex_color_per_cell(*v1, *spheres, cgogn::ui::PointCloudRender<Points>::AttributePerCell::PER_VERTEX);

	NonManifold* skeleton = mpnm.mesh(mps.mesh_name(*s) + "_skeleton");
	auto nm_vertex_position = cgogn::get_attribute<Vec3, NMVertex>(*skeleton, "position");

	srnm.set_vertex_position(*v1, *skeleton, nm_vertex_position);
	srnm.set_edge_width(*v1, *skeleton, 4.0f);

	return app.launch();
}
