/*******************************************************************************
* CGoGN                                                                        *
* Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_MODULE_SHALLOW_WATER_H_
#define CGOGN_MODULE_SHALLOW_WATER_H_

#include <cgogn/ui/module.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class ShallowWater : public ViewModule
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Face = typename mesh_traits<MESH>::Face;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	ShallowWater(const App& app) :
		ViewModule(app, "ShallowWater (" + std::string{mesh_traits<MESH>::name} + ")")
	{}
	~ShallowWater()
	{}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

    void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});


		
		ImGui::End();
	}

private:

	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHALLOW_WATER_H_
