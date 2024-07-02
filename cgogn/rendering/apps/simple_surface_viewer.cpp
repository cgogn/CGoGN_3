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
#include <cgogn/rendering/ui_modules/surface_tex_render.h>




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
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

using Vec3 = cgogn::geometry::Vec3;
using Vec2 = cgogn::geometry::Vec2;
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
	cgogn::ui::SurfaceTexRender<Mesh> str(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&str);


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

		using Vec4 = cgogn::geometry::Vec4;

		// // attrinut de couleur sur l'arête
		// std::shared_ptr<Attribute<Vec3>> edge_color = cgogn::add_attribute<Vec3, Edge>(*m, "edge_color");
		// std::shared_ptr<Attribute<Vec3>> face_color = cgogn::add_attribute<Vec3, Face>(*m, "face_color");


		// plan diagonal de la BB
		// std::pair<Vec3, Vec3> bb = mp.meshes_bb();
		// Vec3 N = (bb.second -  bb.first).normalized();
		// Vec4 plane(N[0],N[1],N[2], - N.dot((bb.first+ bb.second)/Scalar(2)));
		// plane.topRows<3>() = bb.second - bb.first;
	
		// //parcours de arêtes rouge si coupee par le plan, inon bleu 
		// cgogn::foreach_cell(*m, [&](Edge e) -> bool {
		// 	if (cgogn::is_incident_to_boundary(*m,e))
		// 	{
		// 		cgogn::value<Vec3>(*m, edge_color, e) = Vec3(0,1,0);
		// 		std::vector<Face> incF = incident_faces(*m, e);
		// 		for(auto f: incF)
		// 			cgogn::value<Vec3>(*m, face_color, f) = Vec3(1,0,1);
		// 	}
		// 	else
		// 	{
		// 		std::vector<Vertex> iv = incident_vertices(*m, e);
		// 		const Vec3& A = cgogn::value<Vec3>(*m, vertex_position, iv[0]);
		// 		const Vec3& B = cgogn::value<Vec3>(*m, vertex_position, iv[1]);
		// 		bool Sa = std::signbit(plane.dot(Vec4(A[0],A[1],A[2],1)));
		// 		bool Sb = std::signbit(plane.dot(Vec4(B[0],B[1],B[2],1)));
		// 		if (Sa != Sb)
		// 			cgogn::value<Vec3>(*m, edge_color, e) = Vec3(1,0,0);
		// 		else
		// 			cgogn::value<Vec3>(*m, edge_color, e) = Vec3(0,0,1);	
		// 	}
		// 	return true;
		// });


		mp.set_mesh_bb_vertex_position(*m, vertex_position);

		sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());

		sr.set_vertex_position(*v1, *m, vertex_position);
		sr.set_vertex_normal(*v1, *m, vertex_normal);

		str.set_vertex_position(*v1, *m, vertex_position);
		std::shared_ptr<Attribute<Vec2>> vertex_tc = cgogn::add_attribute<Vec2, Vertex>(*m, "tc");	
		
		std::pair<Vec3,Vec3> bb = mp.meshes_bb();
		Vec3 bbw = bb.second - bb.first;

		cgogn::foreach_cell(*m, [&](Vertex v) -> bool {
			Vec2 P = cgogn::value<Vec3>(*m, vertex_position, v).topRows<2>();
			Vec2 TC{(P.x() - bb.first.x())/ bbw.x(), (P.y() - bb.first.y())/ bbw.y()};
			cgogn::value<Vec2>(*m, vertex_tc, v) = TC;
			// std::cout<< TC.transpose() << std::endl;
			return true;
		});
		str.set_vertex_texcoord(*v1, *m, vertex_tc);
		str.chekered_texture();
	}

	return app.launch();
}
