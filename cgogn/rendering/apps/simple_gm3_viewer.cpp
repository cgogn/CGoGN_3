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

#include <cgogn/core/types/cmap/gmap/gmap3.h>
//#include <cgogn/core/types/cmap/cmap3.h>


#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/volume_render.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/functions/distance.h>

#include <cgogn/modeling/algos/volume_utils.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Mesh = cgogn::GMap3;
//using Mesh = cgogn::CMap3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;

int main(int argc, char** argv)
{
	using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
	using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
	using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	std::string filename;
	if (argc >= 2)
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple volume viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::VolumeRender<Mesh> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&vr);

	Mesh* m = nullptr;
	std::shared_ptr<Attribute<Vec3>> vertex_position = nullptr;

	auto setPosV = [&](cgogn::Dart d, const Vec3& P) { cgogn::value<Vec3>(*m, vertex_position, Vertex(d)) = P; };

	if (filename.empty())
	{
		m = mp.add_mesh("simple");
		cgogn::init_cells_indexing<Vertex>(*m);
		vertex_position = cgogn::add_attribute<Vec3, Vertex>(*m, "position");

		cgogn::Dart d_pyra = add_pyramid(*m, 4,false).dart;
		cgogn::Dart d_hexa = add_prism(*m, 4,false).dart;
		for (int i = 0; i < 4; ++i)
		{
			cgogn::phi3_sew(*m, d_pyra, d_hexa);
			d_pyra = cgogn::phi1(*m, d_pyra);
			d_hexa = cgogn::phi_1(*m, d_hexa);
		}


		cgogn::close(*m, false);

		cgogn::foreach_cell(
			*m,
			[&](Vertex v) {
				cgogn::set_index(*m, v, cgogn::new_index<Vertex>(*m));
				return true;
			},
			cgogn::MapBase::TraversalPolicy::DART_MARKING);


		setPosV(d_pyra, Vec3(-1, -1, -1));
		d_pyra = cgogn::phi1(*m, d_pyra);
		setPosV(d_pyra, Vec3(-1, 1, -1));
		d_pyra = cgogn::phi1(*m, d_pyra);
		setPosV(d_pyra, Vec3(1, 1, -1));
		d_pyra = cgogn::phi1(*m, d_pyra);
		setPosV(d_pyra, Vec3(1, -1, -1));
		setPosV(cgogn::phi<2, -1>(*m, d_pyra), Vec3(0, 0, 1));

		cgogn::Dart dh = cgogn::phi<2, 1, 1, 2>(*m, d_hexa);
		setPosV(dh, Vec3(-1, -1, -2));
		dh = cgogn::phi1(*m, dh);
		setPosV(dh, Vec3(-1, 1, -2));
		dh = cgogn::phi1(*m, dh);
		setPosV(dh, Vec3(1, 1, -2));
		dh = cgogn::phi1(*m, dh);
		setPosV(dh, Vec3(1, -1, -2));

//		cgogn::dump_map_darts(*m);

				std::vector<Edge> ve;
		cgogn::foreach_cell(*m, [&](Edge e) {
			ve.push_back(e);
			return true;
		});

		for (Edge e : ve)
		{
			Vertex v1(e.dart);
			Vertex v2(cgogn::phi1(*m,e.dart));
			Vertex v3 = cgogn::cut_edge(*m, e);
			cgogn::value<Vec3>(*m, vertex_position, v3) =
				(cgogn::value<Vec3>(*m, vertex_position, v1) + cgogn::value<Vec3>(*m, vertex_position, v2)) / 2.0;
		}

		dh = d_pyra;
		for (int i=0; i<4; ++i)
		{
			cgogn::Dart dv1 = cgogn::phi<2,1,1>(*m, dh);
			cgogn::Dart dv2 = cgogn::phi<1,1>(*m, dv1);
			cgogn::cut_face(*m, Vertex(dv1), Vertex(dv2));
			dh = cgogn::phi<1, 1>(*m, dh);
		}

		std::vector<cgogn::Dart> vp;
		dh = cgogn::phi<2, 1, 1>(*m, d_pyra);
		vp.push_back(dh);
		dh = cgogn::phi<1, 2, 1>(*m, dh);
		vp.push_back(dh);
		dh = cgogn::phi<1, 2, 1>(*m, dh);
		vp.push_back(dh);
		dh = cgogn::phi<1, 2, 1>(*m, dh);
		vp.push_back(dh);

		cgogn::cut_volume(*m, vp);
	
		cgogn::index_cells<Volume>(*m);

		cgogn::dump_map_darts(*m);

		mp.emit_connectivity_changed(*m);
	}
	else
	{
		m = mp.load_volume_from_file(filename);
		if (!m)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}
		vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	}

	
	vr.set_vertex_position(*v1, *m, vertex_position);
	mp.set_mesh_bb_vertex_position(*m, vertex_position);


	return app.launch();
}
