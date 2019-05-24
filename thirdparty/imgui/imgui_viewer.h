
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

#ifndef CGOGN_IMGUI_VIEWER_H_
#define CGOGN_IMGUI_VIEWER_H_


#include <cgogn/rendering_pureGL/pure_gl_viewer.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace cgogn
{

class ImGUIViewer: public rendering_pgl::PureGLViewer
{
protected:
	bool need_draw_;
	GLFWwindow* window;
	std::string win_name_;
public:
	ImGUIViewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ImGUIViewer);
	ImGUIViewer(int32 w, int32 h);
	~ImGUIViewer();
	virtual void resize_event(int32 w, int32 h);
	virtual void close_event();
	virtual bool init()=0;
	virtual void draw()=0;
	virtual void interface();

	void set_window_title(const std::string&  name);

	bool launch();
};

}
#endif
