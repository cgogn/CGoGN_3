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
#include <chrono>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/modeling/ui_modules/action_unit_change.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/rendering/ui_modules/surface_obj_render.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
// #include <cgogn/rendering/ui_modules/vector_per_vertex_render.h>

#include <cgogn/core/types/mesh_views/cell_filter.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/algos/blending.h>

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

	using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

	using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
	using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
	using Face = typename cgogn::mesh_traits<Mesh>::Face;

	using Vec4 = cgogn::geometry::Vec4;
	using Vec3 = cgogn::geometry::Vec3;
	using Vec2 = cgogn::geometry::Vec2;
	using Scalar = cgogn::geometry::Scalar;

	std::string dirname;
	std::string path_openface;
	if (argc < 5){
		std::cout << "How to launch : Folder with AUs (based on meshes path) , file with texture (based on texture path) , OpenFace path, file with texture normals (based on texture path)\n" << std::endl;
		return 1;
	}
	else{
		dirname = std::string(DEFAULT_MESH_PATH) + std::string(argv[1]);
		path_openface = std::string(argv[3]);
	}
		

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Action Unit Change");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::ActionUnitChange<Mesh> auc(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);
	cgogn::ui::SurfaceModeling<Mesh> sm(app);
	//cgogn::ui::VectorPerVertexRender<Mesh> vpvr(app);
	cgogn::ui::SurfaceObjRender<Mesh> sor(app);

	auc.set_directory(dirname);
	auc.set_pathOpenface(path_openface);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sor);
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
	std::shared_ptr<Attribute<Vec3>> vertex_distance = cgogn::add_attribute<Vec3, Vertex>(*m_pos, "distance");

	auc.set_mesh(*m_pos,vertex_position);
	auc.set_position_attr_name("position");
	auc.set_attribute(*m_pos,vertex_position.get(),"position_interpolation",1.);
	auc.set_attribute(*m_pos,vertex_position.get(),"distance",1.);

	if (argc == 7)
	{
		std::string csv_path_video = argv[5];
		std::string path_video = argv[6];
		auc.exec_mode(true,csv_path_video,path_video);
	}

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
	{
		sor.load_texture(std::string(DEFAULT_TEXTURE_PATH) + std::string(argv[2]));
		sor.load_texture_norm(std::string(DEFAULT_TEXTURE_PATH) + std::string(argv[4]));
	}


	std::vector<std::string> args(argv, argv+argc);
  	for (size_t i = 1; i < args.size(); ++i) {
    	if (args[i] == "-perf") {
			std::ofstream perfFile("performance.txt");
			if (perfFile.is_open())
			{
				int nb_iter = 5000;
				std::vector<std::shared_ptr<Attribute<Vec3>>> attribute_to_blend_;
				std::vector<float> weight_list;

				perfFile << "Performance tests for blending\n\n";

				for (int j = 1; j < auc.pos_aus_.size() ; j++)
				{
					attribute_to_blend_.push_back(auc.pos_aus_[j]);
					
					auto t1 = high_resolution_clock::now();
					for (int k = 0; k < nb_iter; k++)
					{
						for (int z = 0; z < attribute_to_blend_.size(); z++)
						{
							int weight = (std::rand()%5) + 1;
							weight_list.push_back(weight);
						}
						auc.blending(*m_pos,attribute_to_blend_,weight_list);
						weight_list.clear();
					}
					auto t2 = high_resolution_clock::now();

					duration<double, std::milli> ms_double = t2 - t1;

					perfFile << ms_double.count() << "ms for " << nb_iter << " blendings with " << j << " AU with weight change for each AUs\n";
					perfFile << ms_double.count() / float(nb_iter) << "ms per operation for blending with " << j << " AU with weight change for each AUs\n";
				}

				perfFile << "\n\n";

				attribute_to_blend_.clear();
				weight_list.clear();

				for (int j = 1; j < auc.pos_aus_.size() ; j++)
				{
					attribute_to_blend_.push_back(auc.pos_aus_[j]);
					weight_list.push_back((std::rand()%5) + 1);
					
					auto t1 = high_resolution_clock::now();
					for (int k = 0; k < nb_iter; k++)
					{
						auc.blending(*m_pos,attribute_to_blend_,weight_list);
					}
					auto t2 = high_resolution_clock::now();

					duration<double, std::milli> ms_double = t2 - t1;

					perfFile << ms_double.count() << "ms for " << nb_iter << " blendings with " << j << " AU\n";
					perfFile << ms_double.count() / float(nb_iter) << "ms per operation for blending with " << j << " AU\n";

				}

				perfFile << "\n\n";

				attribute_to_blend_.clear();
				weight_list.clear();

				for (int j = 1; j < auc.pos_aus_.size() ; j++)
				{
					attribute_to_blend_.push_back(auc.pos_aus_[j]);
					auto t1 = high_resolution_clock::now();
					for (int k = 0; k < nb_iter; k++)
					{
						for (int z = 0; z < attribute_to_blend_.size(); z++)
						{
							int weight = (std::rand()%5) + 1;
							weight_list.push_back(weight);
						}
						cgogn::modeling::blending(*m_pos,attribute_to_blend_,weight_list,"position");
						weight_list.clear();
					}
					auto t2 = high_resolution_clock::now();

					duration<double, std::milli> ms_double = t2 - t1;

					perfFile << ms_double.count() << "ms for " << nb_iter << " blendings with modeling::blending with " << j << " AU weight change for each AUs\n";
					perfFile << ms_double.count() / float(nb_iter) << "ms per operation for blending with modeling::blending with " << j << " AU weight change for each AUs\n";
				}

				perfFile << "\n\n";

				attribute_to_blend_.clear();
				weight_list.clear();

				for (int j = 1; j < auc.pos_aus_.size() ; j++)
				{
					attribute_to_blend_.push_back(auc.pos_aus_[j]);
					weight_list.push_back(rand()%5 + 1);
					auto t1 = high_resolution_clock::now();
					for (int k = 0; k < nb_iter; k++)
					{
						cgogn::modeling::blending(*m_pos,attribute_to_blend_,weight_list,"position");
					}
					auto t2 = high_resolution_clock::now();

					duration<double, std::milli> ms_double = t2 - t1;

					perfFile << ms_double.count() << "ms for " << nb_iter << " blendings with modeling::blending with " << j << " AU\n";
					perfFile << ms_double.count() / float(nb_iter) << "ms per operation for blending with modeling::blending with " << j << " AU\n";
				}
			}
			perfFile.close();
    	}
  	}
	
	return app.launch();
}
