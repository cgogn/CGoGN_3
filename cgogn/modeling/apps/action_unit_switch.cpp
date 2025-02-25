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

// #define USE_GMAP

#ifdef USE_GMAP
#include <cgogn/core/types/maps/gmap/gmap2.h>
#else
#include <cgogn/core/types/maps/cmap/cmap2.h>
#endif

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui_internal.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/vector_per_vertex_render.h>
#include <cgogn/modeling/ui_modules/action_unit_change.h>
// #include <cgogn/geometry/ui_modules/surface_selection.h>

#include <cgogn/core/types/mesh_views/cell_filter.h>
#include <cgogn/modeling/algos/subdivision.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

#ifdef USE_GMAP
using Mesh = cgogn::GMap2;
#else
using Mesh = cgogn::CMap2;
#endif

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;

int main(int argc, char** argv)
{
	using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
	using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
	using Face = typename cgogn::mesh_traits<Mesh>::Face;

	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	std::string dirname;
	if (argc < 2)
		dirname = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		dirname = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Action Unit Change");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::VectorPerVertexRender<Mesh> vpvr(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);
	cgogn::ui::SurfaceModeling<Mesh> sm(app);
	cgogn::ui::ActionUnitChange<Mesh> auc(app);
	// cgogn::ui::SurfaceSelection<Mesh> ss(app);


	auc.get_directory(dirname);

	Mesh* m = mp.load_surface_from_file(dirname + std::string("AU00.obj"));
	if (!m)
	{
		std::cout << "Folder with files not found" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> vertex_position_interpolation = cgogn::add_attribute<Vec3, Vertex>(*m, "position_interpolation");
	std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");
	std::shared_ptr<Attribute<Vec3>> vertex_color = cgogn::add_attribute<Vec3, Vertex>(*m, "color");
	
	auc.get_mesh(*m,vertex_position);

	

	ImGui::GetIO().DeltaTime = 1./30.;

	cgogn::ui::View* v1 = app.current_view();

	app.init_modules();
	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&vpvr);
	// v1->link_module(&ss);

	mp.set_mesh_bb_vertex_position(*m, vertex_position);

	auc.get_view(*v1);

	sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());

	// mp.emit_connectivity_changed(m);
	// mp.emit_attribute_changed(m, vertex_position.get());

	sr.set_vertex_position(*v1, *m, vertex_position);
	sr.set_vertex_normal(*v1, *m, vertex_normal);

	return app.launch();
}
