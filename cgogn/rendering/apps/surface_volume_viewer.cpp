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

#include <cgogn/core/functions/attributes.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/surface_render_vector.h>
#include <cgogn/rendering/ui_modules/volume_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using SurfaceMesh = cgogn::CMap2;
using VolumeMesh = cgogn::CMap3;

template <typename T>
using SurfaceAttribute = typename cgogn::mesh_traits<SurfaceMesh>::Attribute<T>;
using SurfaceVertex = typename cgogn::mesh_traits<SurfaceMesh>::Vertex;

template <typename T>
using VolumeAttribute = typename cgogn::mesh_traits<VolumeMesh>::Attribute<T>;
using VolumeVertex = typename cgogn::mesh_traits<VolumeMesh>::Vertex;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string surface_filename, volume_filename;
	if (argc < 3)
	{
		volume_filename = std::string(DEFAULT_MESH_PATH) + std::string("tet/hand.tet");
		// std::cout << "Usage: " << argv[0] << " surface_filename volume_filename" << std::endl;
		// return 1;
	}
	else
	{
		surface_filename = std::string(argv[1]);
		volume_filename = std::string(argv[2]);
	}

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<SurfaceMesh> mps(app);
	cgogn::ui::SurfaceRender<SurfaceMesh> sr(app);
	cgogn::ui::SurfaceRenderVector<SurfaceMesh> srv(app);
	cgogn::ui::SurfaceDifferentialProperties<SurfaceMesh> sdp(app);

	cgogn::ui::MeshProvider<VolumeMesh> mpv(app);
	cgogn::ui::VolumeRender<VolumeMesh> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mps);
	v1->link_module(&sr);
	v1->link_module(&srv);

	cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&mpv);
	v2->link_module(&vr);

	SurfaceMesh* sm = mps.load_surface_from_file(surface_filename);
	if (!sm)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<SurfaceAttribute<Vec3>> vertex_position_s =
		cgogn::get_attribute<Vec3, SurfaceVertex>(*sm, "position");
	std::shared_ptr<SurfaceAttribute<Vec3>> vertex_normal_s = cgogn::add_attribute<Vec3, SurfaceVertex>(*sm, "normal");

	sdp.compute_normal(*sm, vertex_position_s.get(), vertex_normal_s.get());

	sr.set_vertex_position(*v1, *sm, vertex_position_s);
	sr.set_vertex_normal(*v1, *sm, vertex_normal_s);

	srv.set_vertex_position(*v1, *sm, vertex_position_s);
	srv.set_vertex_vector(*v1, *sm, vertex_normal_s);

	VolumeMesh* vm = mpv.load_volume_from_file(volume_filename);
	if (!vm)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<VolumeAttribute<Vec3>> vertex_position_v =
		cgogn::get_attribute<Vec3, VolumeVertex>(*vm, "position");

	vr.set_vertex_position(*v2, *vm, vertex_position_v);

	return app.launch();
}
