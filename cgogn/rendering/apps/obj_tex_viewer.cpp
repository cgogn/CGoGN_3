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


#include <cgogn/core/types/maps/cmap/cmap2.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/surface_obj_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Mesh = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

using cgogn::geometry::Vec3;
using cgogn::geometry::Vec2;
using cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string filename;
	std::string imgname;
	if (argc <= 1)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("obj/cube.obj"); 
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Obj tex viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceObjRender<Mesh> sor(app);

	app.init_modules();
	// load texture after init_modules
	if (argc <= 2)
	{
		cgogn::rendering::GLImage img(16, 16, 3);
		std::vector<cgogn::uint8> pix;
		pix.reserve(16*16*3);
		for (int i = 0; i < 16; ++i)
		{
			for (int j = 0; j < 16; ++j)
			{
				if ((i + j) % 2 == 0)
				{
					pix.push_back(i*8);
					pix.push_back(0);
					pix.push_back(j*8);
				}
				else
				{
					pix.push_back(255);
					pix.push_back(255);
					pix.push_back(255);
				}
			}
		}
		img.copy_pixels_data(pix);
		sor.load_texture(img);
	}
	else
		sor.load_texture(std::string(argv[2]));


	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sor);

	if (filename.length() > 0)
	{
		auto [m_pos, m_tc, m_no] = mp.load_surface_from_OBJ_file(filename);
		if (!m_pos)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}
	}

	return app.launch();
}
