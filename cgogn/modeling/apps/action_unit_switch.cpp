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

#ifdef USE_GMAP
#include <cgogn/core/types/maps/gmap/gmap2.h>
#else
#include <cgogn/core/types/maps/cmap/cmap2.h>
#endif

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/io/utils.h>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui_internal.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/modeling/ui_modules/action_unit_change.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/rendering/ui_modules/surface_obj_render.h>
//#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
// #include <cgogn/rendering/ui_modules/vector_per_vertex_render.h>

#include <cgogn/core/types/mesh_views/cell_filter.h>
#include <cgogn/modeling/algos/subdivision.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "meshes/"
#define DEFAULT_TEXTURE_PATH CGOGN_STR(CGOGN_DATA_PATH) "textures/"

using namespace cgogn::numerics;

#ifdef USE_GMAP
using Mesh = cgogn::GMap2;
#else
using Mesh = cgogn::CMap2;
#endif

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;

int main(int argc, char** argv)
{
	using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
	using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
	using Face = typename cgogn::mesh_traits<Mesh>::Face;

	using Vec3 = cgogn::geometry::Vec3;
	using Vec2 = cgogn::geometry::Vec2;
	using Scalar = cgogn::geometry::Scalar;

	std::string dirname;
	if (argc < 2){
		std::cout << "Folder with files not found" << std::endl;
		return 1;
	}
	else
		dirname = std::string(DEFAULT_MESH_PATH) + std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Action Unit Change");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::ActionUnitChange<Mesh> auc(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);
	//cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::SurfaceModeling<Mesh> sm(app);
	//cgogn::ui::VectorPerVertexRender<Mesh> vpvr(app);
	cgogn::ui::SurfaceObjRender<Mesh> sor(app);

	auc.set_directory(dirname);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sor);
	//v1->link_module(&sr);
	v1->link_module(&auc);
	//v1->link_module(&vpvr);

	auto [m_pos,m_tc,m_no] = mp.load_surface_from_OBJ_file(dirname + std::string("AU00.obj"));
	if (!m_pos)
	{
		std::cout << "Folder with files not found" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m_pos, "position");
	std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m_pos, "normal");
	std::shared_ptr<Attribute<Vec3>> vertex_position_interpolation = cgogn::add_attribute<Vec3, Vertex>(*m_pos, "position_interpolation");
	std::shared_ptr<Attribute<Vec3>> vertex_color = cgogn::add_attribute<Vec3, Vertex>(*m_pos, "color");
	std::shared_ptr<Attribute<Vec3>> vertex_distance = cgogn::add_attribute<Vec3, Vertex>(*m_pos, "distance");

	auc.set_mesh(*m_pos,vertex_position);
	auc.set_to_blue(*m_pos);
	auc.set_attribute(*m_pos,vertex_position.get(),"position_interpolation",1.);
	auc.set_attribute(*m_pos,vertex_position.get(),"distance",1.);

	sdp.compute_normal(*m_pos, vertex_position.get(), vertex_normal.get());

	mp.set_mesh_bb_vertex_position(*m_pos, vertex_position);

	auc.set_view(*v1);
	auc.setup_mesh_attributes();
	//auc.setup_csv_matrix();

	if (argc <= 2)
	{
		cgogn::rendering::GLImage img(16, 16, 3);
		std::vector<std::array<cgogn::uint8, 3>> pix;
		pix.reserve(16*16);
		for (int i = 0; i < 16; ++i)
			for (int j = 0; j < 16; ++j)
				if ((i + j) % 2 == 0)
					pix.push_back({0u, 0u, 0u});
				else
					pix.push_back({255u, 255u, 255u});
		img.copy_pixels_data(pix.data()->data());
		sor.load_texture(img);
	}
	else
		sor.load_texture(std::string(DEFAULT_TEXTURE_PATH) + std::string(argv[2]));

	// sr.set_vertex_position(*v1, *m_pos, vertex_position);
	// sr.set_vertex_normal(*v1, *m_pos, vertex_normal);

	
	return app.launch();
}
