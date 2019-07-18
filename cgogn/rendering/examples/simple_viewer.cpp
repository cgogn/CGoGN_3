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

#include <cgogn/core/functions/attributes.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/cmap_provider/cmap_provider.h>
#include <cgogn/ui/modules/surface_differential_properties/surface_differential_properties.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Mesh = cgogn::CMap2;

template <typename T>
using AttributePtr = typename cgogn::mesh_traits<Mesh>::AttributePtr<T>;
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
	app.set_window_title("Simple viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::CMapProvider cmap_provider(app);
	Mesh* m = cmap_provider.import_surface_from_file(filename);

	cgogn::ui::SurfaceRender sr(app);
	cgogn::ui::SurfaceDifferentialProperties sdp(app);

	sr.init();
	sdp.init();

	AttributePtr<Vec3> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	AttributePtr<Vec3> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");
	sdp.compute_normal(*m, vertex_position, vertex_normal);
	sr.update(m, vertex_position, vertex_normal);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&sr);

	return app.launch();
}
