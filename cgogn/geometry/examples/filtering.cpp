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
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_differential_properties/surface_differential_properties.h>
#include <cgogn/ui/modules/surface_filtering/surface_filtering.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>
#include <cgogn/ui/modules/surface_render_vector/surface_render_vector.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

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
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Filtering");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::SurfaceRenderVector<Mesh> srv(app);
	cgogn::ui::SurfaceFiltering<Mesh> sf(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);

	app.init_modules();

	Mesh* m = mp.load_surface_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");
	sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());
	
	sr.set_vertex_position(*m, vertex_position);
	sr.set_vertex_normal(*m, vertex_normal);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&sr);
	v1->link_module(&srv);

	// cgogn::ui::View* v2 = app.add_view();
	// v2->link_module(&sr);
	// v2->link_module(&srv);

	cgogn::ui::MeshData<Mesh>* md = mp.mesh_data(m);
	md->set_bb_vertex_position(vertex_position);
	Vec3 diagonal = md->bb_max_ - md->bb_min_;
	Vec3 center = (md->bb_max_ + md->bb_min_) / 2.0f;
	
	v1->set_scene_radius(diagonal.norm() / 2.0f);
	v1->set_scene_center(center);
	// v2->set_scene_radius(diagonal.norm() / 2.0f);
	// v2->set_scene_center(center);

	return app.launch();
}
