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

#include<cgogn/rendering/shaders/shader_function_color_maps.h>

#ifndef CGOGN_MODULE_TOOLS_H_
#define CGOGN_MODULE_TOOLS_H_

#include <imgui/imgui.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/rendering/shaders/shader_function_color_maps.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/attributes.h>

namespace cgogn
{
namespace ui
{

/**
 * @brief generate combo for attribute selection
 * @param m mesh
 * @param att attribute
 * @param f code to execute when combo selection change
 */
template <typename CELL, typename T, typename MESH, typename FUNC>
inline void imgui_combo_attribute(const MESH& m,
					std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& att,
					const std::string& label,
					const FUNC& f)
{
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	std::string selected_attrib = att ? att->name() : "-- select --";
	if (ImGui::BeginCombo(label.c_str(), selected_attrib.c_str()))
	{
		foreach_attribute<T,CELL>(m,
				  [&](const std::shared_ptr<Attribute>& attribute)
		{
			bool is_selected = attribute == att;
			if (ImGui::Selectable(attribute->name().c_str(), is_selected))
			{
				f(attribute);
			}
			if (is_selected)
				ImGui::SetItemDefaultFocus();
		});
		ImGui::EndCombo();
	}

	if (att)
	{
		double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
		ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
		if (ImGui::Button((std::string("X##")+label).c_str()))
		{
			f(nullptr);
		}
	}
}

//if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
//{
//	foreach_attribute<Vec3, Vertex>(
//		*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
//			bool is_selected = attribute == p.vertex_position_;
//			if (ImGui::Selectable(attribute->name().c_str(), is_selected))
//				set_vertex_position(*selected_view_, *selected_mesh_, attribute);
//			if (is_selected)
//				ImGui::SetItemDefaultFocus();
//		});
//	ImGui::EndCombo();
//}
//if (p.vertex_position_)
//{
//	ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
//	if (ImGui::Button("X##position"))
//		set_vertex_position(*selected_view_, *selected_mesh_, nullptr);
//}





bool imgui_colormap_interface(rendering::shader_funcion::ColorMap::Uniforms& cm, const std::string& label)
{
	bool need_update = false;
	ImGui::TextUnformatted("ColorMAP:");
	ImGui::BeginGroup();
	need_update |= ImGui::RadioButton((std::string("BWR##")+label).c_str(), &(cm.color_map_),0);ImGui::SameLine();
	need_update |= ImGui::RadioButton((std::string("CWR##")+label).c_str(), &cm.color_map_, 1);ImGui::SameLine();
	need_update |= ImGui::RadioButton((std::string("BCGYR##")+label).c_str(), &cm.color_map_, 2);ImGui::SameLine();
	need_update |= ImGui::RadioButton((std::string("BCR##")+label).c_str(), &cm.color_map_, 3);
	need_update |= ImGui::SliderInt((std::string("expansion##")+label).c_str(), &(cm.expansion_), -2, 2);
	need_update |= ImGui::InputFloat((std::string("min##")+label).c_str(), &(cm.min_value_));
	need_update |= ImGui::InputFloat((std::string("max##")+label).c_str(), &(cm.max_value_));
	ImGui::EndGroup();
	return need_update;
}

}
}

#endif // CGOGN_MODULE_MESH_PROVIDER_H_
