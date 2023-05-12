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

#define USE_GMAP

#ifdef USE_GMAP
#include <cgogn/core/types/maps/gmap/gmap2.h>
#else
#include <cgogn/core/types/maps/cmap/cmap2.h>
#endif


#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/modeling/algos/subdivision/surface_catmull_clark.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/core/types/mesh_views/cell_filter.h>
#include <cgogn/core/types/mesh_traits.h>

//#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

#ifdef USE_GMAP
using Mesh = cgogn::GMap2;
#else
using Mesh = cgogn::CMap2;
#endif

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;
using namespace cgogn::numerics;

Volume ring(Mesh& m, uint32 nb, double r, double R)
{
	using PMesh = typename cgogn::mesh_traits<Mesh>::Parent;
	using Face = typename cgogn::mesh_traits<Mesh>::Face;

	auto e0 = cgogn::add_face(static_cast<PMesh&>(m), 4, false);
	cgogn::Dart dp = cgogn::phi<1, 1>(m, e0.dart);
	for (uint32 i = 1; i < nb; ++i)
	{
		auto e = cgogn::add_face(static_cast<PMesh&>(m), 4, false);
		cgogn::phi2_sew(m, dp, e.dart);
		dp = cgogn::phi<1, 1>(m, e.dart);
	}
	cgogn::phi2_sew(m, dp, e0.dart);

	Face h = cgogn::close_hole(m, cgogn::phi_1(m, e0.dart), true);
	cgogn::foreach_dart_of_orbit(m, h, [&](cgogn::Dart hd) -> bool {
		set_boundary(m, hd, true);
		return true;
	});

	h = cgogn::close_hole(m, cgogn::phi1(m, e0.dart), true);
	cgogn::foreach_dart_of_orbit(m, h, [&](cgogn::Dart hd) -> bool {
		set_boundary(m, hd, true);
		return true;
	});

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");

	cgogn::Dart d = e0.dart;
	for (uint32 i = 0; i < nb; ++i)
	{
		double alpha = i * (2.0 * M_PI / nb);

		Vertex v(d);
		cgogn::set_index<Vertex>(m, v, cgogn::new_index<Vertex>(m));
		cgogn::value<Vec3>(m, vertex_position, v) = Vec3(r * std::cos(alpha), r * std::sin(alpha), -0.5);
		d = phi1(m, d);
		v = Vertex(d);
		cgogn::set_index<Vertex>(m, v, cgogn::new_index<Vertex>(m));
		cgogn::value<Vec3>(m, vertex_position, v) = Vec3(R * std::cos(alpha), R * std::sin(alpha), 0.5);
		d = cgogn::phi<1, 2>(m, d);
	}

	return Volume(phi_1(m, e0.dart));
}

#ifdef USE_GMAP
Volume moebius(Mesh& m, uint32 nb, double r, double R)
{
	using PMesh = typename cgogn::mesh_traits<Mesh>::Parent;
	using Face = typename cgogn::mesh_traits<Mesh>::Face;

	auto e0 = cgogn::add_face(static_cast<PMesh&>(m), 4, false);
	cgogn::Dart dp = cgogn::phi<1, 1>(m, e0.dart);
	for (uint32 i = 1; i < nb; ++i)
	{
		auto e = cgogn::add_face(static_cast<PMesh&>(m), 4, false);
		cgogn::phi2_sew(m, dp, e.dart);
		dp = cgogn::phi<1, 1>(m, e.dart);
	}
	cgogn::phi2_sew(m, dp, beta0(m, e0.dart));

	Face h = cgogn::close_hole(m, cgogn::phi_1(m, e0.dart), false);
	cgogn::foreach_dart_of_orbit(m, h, [&](cgogn::Dart hd) -> bool {
		set_boundary(m, hd, true);
		return true;
	});

	// auto vertex_position = cgogn::get_attribute<Vec3, Vertex>(*mesh_, "position");
	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");

	cgogn::Dart d = e0.dart;
	for (uint32 i = 0; i < nb; ++i)
	{
		double alpha = i * (2.0 * M_PI / nb);

		double a = (r + R) / 2;
		double b = (R - r) / 2;
		double beta = alpha / 2;
		double rad1 = (a + b * std::cos(beta + M_PI));
		double rad2 = (a + b * std::cos(beta));

		Vertex v(d);
		cgogn::set_index<Vertex>(m, v, cgogn::new_index<Vertex>(m));
		cgogn::value<Vec3>(m, vertex_position, v) =
			Vec3(rad1 * std::cos(alpha), rad1 * std::sin(alpha), b * std::sin(beta + M_PI));
		d = phi1(m, d);
		v = Vertex(d);
		cgogn::set_index<Vertex>(m, v, cgogn::new_index<Vertex>(m));
		cgogn::value<Vec3>(m, vertex_position, v) =
			Vec3(rad2 * std::cos(alpha), rad2 * std::sin(alpha), b * std::sin(beta));
		d = cgogn::phi<1, 2>(m, d);
	}

	return Volume(phi_1(m, e0.dart));
}
#endif


class LocalInterface : public cgogn::ui::Module
{

public:
	LocalInterface(const cgogn::ui::App& app )
		: cgogn::ui::Module(app, "LocalInterface"), mesh_(nullptr), vertex_position_(nullptr), mesh_provider_(nullptr),
		  surf_render_(nullptr)
	{
	}

	~LocalInterface()
	{
	}

	void create()
	{
		mesh_ = mesh_provider_->add_mesh("pyrad_and_prism");
		vertex_position_ = cgogn::add_attribute<Vec3, Vertex>(*mesh_, "position");
		vertex_normal_ = cgogn::add_attribute<Vec3, Vertex>(*mesh_, "normal");

		
		cgogn::init_cells_indexing<Vertex>(*mesh_);
		cgogn::Dart d;
		Vec3 center;

		auto setPosV = [&](const Vec3& P) { cgogn::value<Vec3>(*mesh_, vertex_position_, Vertex(d)) = P + center; };

		d_pyr_ = add_pyramid(*mesh_, 4);
		d = d_pyr_.dart;	
		center = Vec3(0, 3, 0);
		setPosV(Vec3(-1, -1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(-1, 1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(1, 1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(1, -1, -1));
		d = cgogn::phi<2, -1>(*mesh_, d);
		setPosV(Vec3(0, 0, 1));

		d_tet_ = add_pyramid(*mesh_, 3);
		d = d_tet_.dart;
		center = Vec3(3, 0, 0);
		setPosV(Vec3(-1.22456, -0.707, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(0, 1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(1.22456, -0.707, -1));
		d = cgogn::phi<2, -1>(*mesh_, d);
		setPosV(Vec3(0, 0, 1));

		d_prism_ = add_prism(*mesh_, 3);
		d = d_prism_.dart;
		center = Vec3(0, -3, 0);
		setPosV(Vec3(-1.22456, -0.707, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(0, 1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(1.22456, -0.707, -1));
		d = cgogn::phi<2, 1, 1, 2>(*mesh_, d);
		setPosV(Vec3(-1.22456, -0.707, 1));
		d = cgogn::phi_1(*mesh_, d);
		setPosV(Vec3(0, 1, 1));
		d = cgogn::phi_1(*mesh_, d);
		setPosV(Vec3(1.22456, -0.707, 1));

		d_hex_ = add_prism(*mesh_, 4);
		d = d_hex_.dart;
		center = Vec3(-3, 0, 0);
		setPosV(Vec3(-1, -1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(-1, 1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(1, 1, -1));
		d = cgogn::phi1(*mesh_, d);
		setPosV(Vec3(1, -1, -1));
		d = cgogn::phi<2, 1, 1, 2>(*mesh_, d);
		setPosV(Vec3(-1, -1, 1));
		d = cgogn::phi_1(*mesh_, d);
		setPosV(Vec3(-1, 1, 1));
		d = cgogn::phi_1(*mesh_, d);
		setPosV(Vec3(1, 1, 1));
		d = cgogn::phi_1(*mesh_, d);
		setPosV(Vec3(1, -1, 1));

		#ifdef USE_GMAP
			d_moeb_ = moebius(*mesh_, 8, 0.7, 1.3);
		#else
			d_moeb_ = ring(*mesh_, 8, 0.7, 1.3);	
		#endif

		//cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position_.get());
		//cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position_.get());
		// cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position.get());

		cgogn::geometry::compute_normal<Vertex>(*mesh_, vertex_position_.get(), vertex_normal_.get());

		surf_render_->set_vertex_position(*app_.current_view(), *mesh_, vertex_position_);
		surf_render_->set_vertex_normal(*app_.current_view(), *mesh_, vertex_normal_);
		mesh_provider_->set_mesh_bb_vertex_position(*mesh_, vertex_position_);

		mesh_provider_->emit_connectivity_changed(*mesh_);
	}


protected:
	void init() override
	{
		mesh_provider_ = static_cast<cgogn::ui::MeshProvider<Mesh>*>(
			app_.module("MeshProvider (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));

		surf_render_ = static_cast<cgogn::ui::SurfaceRender<Mesh>*>(
			app_.module("SurfaceRender (" + std::string{cgogn::mesh_traits<Mesh>::name} + ")"));
	}

	void left_panel() override
	{
		if (ImGui::Button("CatmullClark"))
		{
			cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position_.get());
			cgogn::geometry::compute_normal<Vertex>(*mesh_, vertex_position_.get(), vertex_normal_.get());

			mesh_provider_->emit_connectivity_changed(*mesh_);
			mesh_provider_->emit_attribute_changed(*mesh_, vertex_position_.get());
		}

		if (ImGui::Button("CatmullClark Meobius"))
		{
			cgogn::DartMarker<Mesh> dm(*mesh_);
			cgogn::foreach_dart_of_orbit(*mesh_, Volume(d_moeb_), [&](cgogn::Dart d) {
				dm.mark(d);
				return true;
			});
		
			cgogn::CellFilter<Mesh> cf(*mesh_);
			cf.set_filter<Vertex>([&](Vertex v) { return dm.is_marked(v.dart); });
			cf.set_filter<Edge>([&](Edge v) { return dm.is_marked(v.dart); });
			cf.set_filter<Face>([&](Face v) { 
				return dm.is_marked(v.dart); 
				});

			std::cout << std::boolalpha << std::is_convertible_v<cgogn::CellFilter<Mesh>&, cgogn::MapBase&> << " "
					  << int(cgogn::mesh_traits<cgogn::CellFilter<Mesh>>::dimension) << std::endl;

			cgogn::modeling::subdivide_catmull_clark(cf, vertex_position_.get());
			mesh_provider_->emit_connectivity_changed(*mesh_);
			mesh_provider_->emit_attribute_changed(*mesh_, vertex_position_.get());
		}

	}

private:
	Mesh* mesh_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	std::shared_ptr<Attribute<Vec3>> vertex_normal_;
	cgogn::ui::MeshProvider<Mesh>* mesh_provider_;
	cgogn::ui::SurfaceRender<Mesh>* surf_render_;
	Volume d_moeb_;
	Volume d_pyr_;
	Volume d_tet_;
	Volume d_hex_;
	Volume d_prism_;


};






int main(int argc, char** argv)
{
	std::string filename = "";
	if (argc > 1)
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Simple surface viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	LocalInterface interf(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);

	interf.create();

//	if (filename.length() > 0)
//	{
//		Mesh* m = mp.load_surface_from_file(filename);
//		if (!m)
//		{
//			std::cout << "File could not be loaded" << std::endl;
//			return 1;
//		}
//
//		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*mesh_, "position");
//		std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*mesh_, "normal");
//
//		 std::shared_ptr<Attribute<Vec3>> face_color = cgogn::add_attribute<Vec3, Face>(*mesh_, "color");
//		 std::shared_ptr<Attribute<Scalar>> face_weight = cgogn::add_attribute<Scalar, Face>(*mesh_, "weight");
//
//		 cgogn::foreach_cell(*mesh_, [&](Face f) -> bool {
//		 	Vec3 c(0, 0, 0);
//		 	c[rand() % 3] = 1;
//		 	cgogn::value<Vec3>(*mesh_, face_color, f) = c;
//		 	cgogn::value<Scalar>(*mesh_, face_weight, f) = double(rand()) / RAND_MAX;
//		 	return true;
//		 });
//
//		 mp.set_mesh_bb_vertex_position(*mesh_, vertex_position);
//
////		sdp.compute_normal(*mesh_, vertex_position.get(), vertex_normal.get());
//
//		sr.set_vertex_position(*v1, *mesh_, vertex_position);
////		sr.set_vertex_normal(*v1, *mesh_, vertex_normal);
//	}
//	else
//	{

		//Mesh* m = mp.add_mesh("pyradmid");
		//std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::add_attribute<Vec3, Vertex>(*mesh_, "position");
		//cgogn::init_cells_indexing<Vertex>(*mesh_);

		////auto pyr = add_pyramid(*mesh_, 4);
		////auto tet = add_pyramid(*mesh_, 3);
		////auto hex = add_prism(*mesh_, 4);

		//cgogn::Dart d;
		//Vec3 center;

		//auto setPosV = [&](const Vec3& P) { cgogn::value<Vec3>(*mesh_, vertex_position, Vertex(d)) = P + center; };

		//d = add_pyramid(*mesh_, 4).dart;
		//center = Vec3(0, 3, 0);
		//setPosV(Vec3(-1, -1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(-1, 1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(1, 1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(1, -1, -1));
		//d = cgogn::phi<2, -1>(*mesh_, d);
		//setPosV(Vec3(0, 0, 1));

		//d = add_pyramid(*mesh_, 3).dart;
		//center = Vec3(3, 0, 0);
		//setPosV(Vec3(-1.22456, -0.707, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(0, 1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(1.22456, -0.707, -1));
		//d = cgogn::phi<2, -1>(*mesh_, d);
		//setPosV(Vec3(0, 0, 1));

		//d = add_prism(*mesh_, 3).dart;
		//center = Vec3(0,- 3, 0);
		//setPosV(Vec3(-1.22456, -0.707, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(0, 1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(1.22456, -0.707, -1));
		//d = cgogn::phi<2, 1, 1, 2>(*mesh_, d);
		//setPosV(Vec3(-1.22456, -0.707, 1));
		//d = cgogn::phi_1(*mesh_, d);
		//setPosV(Vec3(0, 1, 1));
		//d = cgogn::phi_1(*mesh_, d);
		//setPosV(Vec3(1.22456, -0.707, 1));


		//d = add_prism(*mesh_, 4).dart;
		//center = Vec3(-3, 0, 0);
		//setPosV(Vec3(-1, -1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(-1, 1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(1, 1, -1));
		//d = cgogn::phi1(*mesh_, d);
		//setPosV(Vec3(1, -1, -1));
		//d = cgogn::phi<2, 1, 1, 2>(*mesh_, d);
		//setPosV(Vec3(-1, -1, 1));
		//d = cgogn::phi_1(*mesh_, d);
		//setPosV(Vec3(-1, 1, 1));
		//d = cgogn::phi_1(*mesh_, d);
		//setPosV(Vec3(1, 1, 1));
		//d = cgogn::phi_1(*mesh_, d);
		//setPosV(Vec3(1, -1, 1));

		//cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position.get());
		//cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position.get());
		//cgogn::modeling::subdivide_catmull_clark(*mesh_, vertex_position.get());

		//moebius(*mesh_, 32, 0.7, 1.3);

		//mp.set_mesh_bb_vertex_position(*mesh_, vertex_position);
		//sr.set_vertex_position(*v1, *mesh_, vertex_position);
		//mp.emit_connectivity_changed(*mesh_);
//	}

	return app.launch();
}
