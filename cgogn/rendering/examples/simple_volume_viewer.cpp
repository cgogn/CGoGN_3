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

#include <cgogn/core/functions/attributes.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/volume_render/volume_render.h>
#include <cgogn/ui/modules/topo_render/topo_render.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/functions/distance.h>

#include <typeinfo>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH)"/meshes/"

using Mesh = cgogn::CMap3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;


using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{

	std::string filename;
	if (argc < 2)
	{
		filename = std::string(DEFAULT_MESH_PATH) + std::string("tet/hex_dominant.meshb");
//		std::cout << "Usage: " << argv[0] << " volume_mesh_file" << std::endl;
//		return 1;
	}
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::Volume_Render<Mesh> vr(app);
	cgogn::ui::Topo_Render<Mesh> tr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
//	cgogn::ui::View* v2 = app.add_view();
//	cgogn::ui::View* v3 = app.add_view();
//	cgogn::ui::View* v4 = app.add_view();

	v1->link_module(&mp);
	v1->link_module(&vr);
	v1->link_module(&tr);

//	v2->link_module(&mp);
//	v2->link_module(&vr);
//	v2->link_module(&tr);
//	v3->link_module(&mp);
//	v3->link_module(&vr);
//	v3->link_module(&tr);
//	v4->link_module(&mp);
//	v4->link_module(&vr);
//	v4->link_module(&tr);


	Mesh* m = mp.load_volume_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> volume_center = cgogn::add_attribute<Vec3, Volume>(*m, "center");

//	cgogn::index_cells<Volume>(*m);
//	cgogn::geometry::compute_centroid<Vec3,Volume>(*m,vertex_position.get(),volume_center.get());


	mp.set_mesh_bb_vertex_position(m, vertex_position);

	vr.set_vertex_position(*v1, *m, vertex_position);
	tr.set_vertex_position(*v1, *m, vertex_position);


	return app.launch();
}
