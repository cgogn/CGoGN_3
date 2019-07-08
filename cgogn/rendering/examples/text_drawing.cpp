/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/rendering/imgui_viewer.h>
#include <cgogn/rendering/text_drawer.h>
#include <cgogn/core/utils/numerics.h>
#include <chrono>
#include <iomanip>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Vec3 = Eigen::Vector3d;

using namespace cgogn;
namespace RGL = rendering;
class App;
class TextDrawing : public RGL::ImGUIViewer
{
	friend class App;
public:
	TextDrawing();
//	TextDrawing(TextDrawing* ptr);
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(TextDrawing);

	void draw() override;
	void init() override;
	void key_press_event(int k) override;
	void resize_event(int32 w, int32 h) override;

	std::unique_ptr<RGL::TextDrawer> tdr_;
	std::unique_ptr<RGL::TextDrawer::Renderer> tdr_rend_;

	std::unique_ptr<RGL::TextDrawer> tdr2_;
	std::unique_ptr<RGL::TextDrawer::Renderer> tdr_rend2_;

	std::chrono::time_point<std::chrono::system_clock> start_fps_;
	int nb_fps_;
	std::string fps_;
};

class App: public RGL::ImGUIApp
{
public:
	App() {}
};

TextDrawing::TextDrawing() :
	tdr_(nullptr),
	tdr_rend_(nullptr),
	fps_("?????")
{}


void TextDrawing::draw()
{
	RGL::GLMat4 proj = get_projection_matrix();
	RGL::GLMat4 view = get_modelview_matrix();

	tdr_rend_->draw(proj, view);
	RGL::GLMat4 Id = RGL::GLMat4::Identity();
	RGL::Transfo3d t= Eigen::Translation<double,3>(-1,-1,0)*Eigen::Scaling(0.5, 0.5*width()/height(),0.0);
	RGL::GLMat4 ratio = t.matrix().cast<float>();
//	ratio.translate(-1,-1,0);
//	ratio.scale(0.5f, 0.5f*width()/height(),0.0f);
	tdr_rend2_->draw(ratio,Id);

	nb_fps_++;
	if (nb_fps_ == 50)
	{
		std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
		std::chrono::duration<float> elapsed_seconds = end - start_fps_;
		float fps = 50.0f/elapsed_seconds.count();
		fps_ = std::to_string(fps).substr(0,5);
		tdr2_->update_text(0,fps_);
		start_fps_ = std::chrono::system_clock::now();
		nb_fps_ = 0;
	}
}


void TextDrawing::init()
{
	set_scene_radius(5.0);
	set_scene_center(Vec3(0,0,0));
	glClearColor(0.1f,0.1f,0.2f,0.0f);

	tdr_ = cgogn::make_unique<RGL::TextDrawer>();
	tdr_rend_ = tdr_->generate_renderer();

	for (float z=-4; z<4; z += 1)
		for (float y = -4; y<4; y += 1)
			for (float x = -4; x < 4; x += 1)
			{
				Vec3 P{ x,y,z };
				RGL::GLColor col(rand()%255, rand() % 255, rand() % 255, 255);
				col/=255.0;
				float sz = 0.1f*rand() / RAND_MAX + 0.05f;
				std::stringstream ss;
				ss << std::setprecision(2)<< "(" << x << "," << y << "," << z << ")";
				*tdr_ << P << col << sz << ss.str();
			}
	*tdr_ << RGL::TextDrawer::end;

	tdr2_ = cgogn::make_unique<RGL::TextDrawer>();
	tdr_rend2_ = tdr2_->generate_renderer();

	float sz = 32.0f / (/*device_pixel_ratio()* */width());
	*tdr2_ <<Vec3{sz,sz,-1} << RGL::GLColor(1,1,1,1) <<sz << fps_ << " fps"<< RGL::TextDrawer::end ;
	tdr_rend2_->set_italic(0.2f);
	start_fps_ = std::chrono::system_clock::now();
	nb_fps_ = 0;
}

void TextDrawing::resize_event(int32 w, int32 h)
{
	if (tdr2_)
	{
		float sz = 32.0f / (/*device_pixel_ratio()* */w);
		*tdr2_ <<Vec3{sz,sz,-1} << RGL::GLColor(1,1,1,1) <<sz << fps_ << " fps"<< RGL::TextDrawer::end ;
	}

	ImGUIViewer::resize_event(w,h);
}

void TextDrawing::key_press_event(int k)
{
	switch(k)
	{
		case int('A'):
				tdr_->update_text(0,"XXXXXXXXXX");
			break;
		case int('-'):
				tdr_->scale_text(0.9f);
			break;
		case int('+'):
				tdr_->scale_text(1.1f);
			break;
		default:
			break;
	}
}


int main(int argc, char** argv)
{
//	qoglviewer::init_ogl_context();
	App app;
	gl3wInit();
	TextDrawing view;
	app.add_view(&view);
	app.set_window_title("TextDrawing");
	app.launch();
}
