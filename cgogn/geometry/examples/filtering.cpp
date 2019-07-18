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

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;

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

	cgogn::ui::MeshProvider<Mesh> mesh_provider(app);
	mesh_provider.import_surface_from_file(filename);

	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::SurfaceFiltering<Mesh> sf(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);

	sr.init();
	sf.init();
	sdp.init();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&sr);

	cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&sr);

	// cgogn::ui::View* v3 = app.add_view();
	// v3->link_module(&f);

	// cgogn::ui::View* v4 = app.add_view();
	// v4->link_module(&f);

	return app.launch();
}
