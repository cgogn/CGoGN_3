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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/types/attribute_handler.h>
//#include <cgogn/core/functions/attributes.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>

//#include <cgogn/ui/modules/topo_render/topo_render.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/functions/distance.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Mesh = cgogn::CMap3;
template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;

int main(int argc, char** argv)
{

	using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
	using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	std::string filename;
	if (argc < 2)
	{
		filename = std::string(DEFAULT_MESH_PATH) + std::string("tet/hex_dominant.meshb");
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple volume viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::Volume_Render<Mesh> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mp);
	v1->link_module(&vr);

	//	cgogn::ui::View* v2 = app.add_view();
	//	v2->link_module(&mp);
	//	v2->link_module(&vr);

	Mesh* m = mp.load_volume_from_file(filename);

	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> volume_center = cgogn::add_attribute<Vec3, Volume>(*m, "center");
	std::shared_ptr<Attribute<Scalar>> volume_scal = cgogn::add_attribute<Scalar, Volume>(*m, "scal");
	std::shared_ptr<Attribute<Vec3>> volume_color = cgogn::add_attribute<Vec3, Volume>(*m, "color");

	// example of convenient functions for nice old-cgogn2 syntax
	//	auto color_handler = [m,&volume_color] (Volume v) -> Vec3&
	//			{ return cgogn::value<Vec3>(*m, volume_color, v); };
	//	auto scalar_handler = [m,&volume_scal] (Volume v) -> Scalar&
	//			{ return cgogn::value<Scalar>(*m, volume_scal, v); };

	auto color_handler = cgogn::attribute_handler<Volume>(m, volume_color);
	auto scalar_handler = cgogn::attribute_handler<Volume>(m, volume_scal);

	cgogn::index_cells<Volume>(*m);
	cgogn::foreach_cell(*m, [&](Volume v) -> bool {
		Vec3 c(0, 0, 0);
		c[rand() % 3] = 1;
		color_handler[v] = c;
		scalar_handler[v] = double(rand()) / RAND_MAX;
		return true;
	});

	cgogn::geometry::compute_centroid<Vec3, Volume>(*m, vertex_position.get(), volume_center.get());

	mp.set_mesh_bb_vertex_position(m, vertex_position);
	//	vr.set_vertex_position(*v1, *m, vertex_position);
	//	vr.set_volume_scalar(*v1, *m, volume_scal);
	//	vr.set_volume_color(*v1, *m, volume_color);

	return app.launch();
}
