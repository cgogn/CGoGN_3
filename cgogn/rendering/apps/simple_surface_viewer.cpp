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

// #ifdef USE_GMAP
// #include <cgogn/core/types/maps/gmap/gmap2.h>
// #else
#include <cgogn/core/types/maps/cmap/cmap2.h>
// #endif

// #include <cgogn/core/types/triangle_soup/triangle_soup.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

// #ifdef USE_GMAP
// using Mesh = cgogn::GMap2;
// #else
using Mesh = cgogn::CMap2;
// #endif
// using Mesh = cgogn::IncidenceGraph;

// using Mesh = cgogn::TriangleSoup;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string filename;
	if (argc <= 1)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple surface viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);

	// cgogn::ui::View* v2 = app.add_view();
	// v2->link_module(&mp);
	// v2->link_module(&sr);

	if (filename.length() > 0)
	{
		Mesh* m = mp.load_surface_from_file(filename);
		if (!m)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}

		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
		std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");

		std::shared_ptr<Attribute<Vec3>> face_color = cgogn::add_attribute<Vec3, Face>(*m, "color");
		std::shared_ptr<Attribute<Scalar>> face_weight = cgogn::add_attribute<Scalar, Face>(*m, "weight");

		cgogn::foreach_cell(*m, [&](Face f) -> bool {
			Vec3 c(0, 0, 0);
			c[rand() % 3] = 1;
			cgogn::value<Vec3>(*m, face_color, f) = c;
			cgogn::value<Scalar>(*m, face_weight, f) = double(rand()) / RAND_MAX;
			return true;
		});

		mp.set_mesh_bb_vertex_position(*m, vertex_position);

		sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());

		sr.set_vertex_position(*v1, *m, vertex_position);
		sr.set_vertex_normal(*v1, *m, vertex_normal);
	}

	return app.launch();
}
