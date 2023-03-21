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

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/volume_render.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/functions/distance.h>

#include <cgogn/modeling/algos/volume_utils.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/io/utils.h>
#include <cgogn/io/volume/volume_import.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<uint32_t, K>;
using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
using Point = Delaunay::Point;

using Mesh = cgogn::CMap3;
template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using namespace cgogn::numerics;

using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

void test_delaunay(Mesh* m)
{
	cgogn::io::VolumeImportData volume_data;

	volume_data.vertex_position_.emplace_back(0, 0, 0);
	volume_data.vertex_position_.emplace_back(1, 0, 0);
	volume_data.vertex_position_.emplace_back(0, 1, 0);
	volume_data.vertex_position_.emplace_back(0, 0, 1);
	volume_data.vertex_position_.emplace_back(1, 1, 1);
	volume_data.vertex_position_.emplace_back(-1, 0, 1);

	uint32 nb_vertices = volume_data.vertex_position_.size();

	std::vector<unsigned> indices;
	indices.reserve(nb_vertices);
	for (int i = 0; i < nb_vertices; ++i)
		indices.push_back(i);

	// risky !!
	std::vector<Point>* points = reinterpret_cast<std::vector<Point>*>(&volume_data.vertex_position_);

	Delaunay tri(boost::make_zip_iterator(boost::make_tuple(points->begin(), indices.begin())),
				 boost::make_zip_iterator(boost::make_tuple(points->end(), indices.end())));

	uint32 nb_volumes = tri.number_of_finite_cells();
	volume_data.reserve(nb_vertices, nb_volumes);

	for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
	{
		volume_data.volumes_types_.push_back(cgogn::io::VolumeType::Tetra);
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(0)->info());
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(1)->info());
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(2)->info());
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(3)->info());
	}

	import_volume_data(*m, volume_data);

	// Why ?
	cgogn::index_cells<Volume>(*m);
}

int main(int argc, char** argv)
{

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Delaunay volume viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::VolumeRender<Mesh> vr(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&vr);

	Mesh* m = mp.add_mesh("delaunay");
	test_delaunay(m);

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	mp.set_mesh_bb_vertex_position(*m, vertex_position);
	vr.set_vertex_position(*v1, *m, vertex_position);

	return app.launch();
}