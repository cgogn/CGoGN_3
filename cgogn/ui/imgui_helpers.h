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

#ifndef CGOGN_UI_IMGUI_HELPERS_H_
#define CGOGN_UI_IMGUI_HELPERS_H_

#include <imgui/imgui.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/rendering/shaders/shader_function_color_maps.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

namespace cgogn
{

namespace ui
{

/**
 * @brief generate combo box for attribute selection
 * @param m mesh
 * @param label label of the combo box
 * @param selected_attribute current selected attribute
 * @param on_change function to call with newly selected attribute
 */
template <typename CELL, typename T, typename MESH, typename FUNC>
void imgui_combo_attribute(const MESH& m,
						   const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& selected_attribute,
						   const std::string& label, const FUNC& on_change)
{
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<Attribute>&>::value,
				  "Wrong function attribute parameter type");

	if (ImGui::BeginCombo(label.c_str(), selected_attribute ? selected_attribute->name().c_str() : "-- select --"))
	{
		foreach_attribute<T, CELL>(m, [&](const std::shared_ptr<Attribute>& attribute) {
			bool is_selected = attribute == selected_attribute;
			if (ImGui::Selectable(attribute->name().c_str(), is_selected))
				on_change(attribute);
			if (is_selected)
				ImGui::SetItemDefaultFocus();
		});
		ImGui::EndCombo();
	}
	if (selected_attribute)
	{
		double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
		ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - float32(X_button_width));
		if (ImGui::Button(("X##" + label).c_str()))
			on_change(nullptr);
	}
}

template <template <typename MESH> typename MESH_PROVIDER, typename MESH, typename FUNC>
bool imgui_mesh_selector(MESH_PROVIDER<MESH>* mesh_provider, const MESH* selected_mesh, const std::string& label,
						 const FUNC& on_change)
{
	static_assert(is_func_parameter_same<FUNC, MESH&>::value, "Wrong function parameter type");
	if (ImGui::ListBoxHeader(label.c_str(), mesh_provider->number_of_meshes()))
	{
		mesh_provider->foreach_mesh([&](MESH& m, const std::string& name) {
			if (ImGui::Selectable(name.c_str(), &m == selected_mesh))
				on_change(m);
		});
		ImGui::ListBoxFooter();
		return true;
	}
	return false;
}

template <typename FUNC>
bool imgui_view_selector(ViewModule* vm, const View* selected, const FUNC& on_change)
{
	static_assert(is_func_parameter_same<FUNC, View*>::value, "Wrong function parameter type");
	if (ImGui::BeginCombo("View", selected->name().c_str()))
	{
		for (View* v : vm->linked_views())
		{
			bool is_selected = v == selected;
			if (ImGui::Selectable(v->name().c_str(), is_selected))
				on_change(v);
			if (is_selected)
				ImGui::SetItemDefaultFocus();
		}
		ImGui::EndCombo();
		return true;
	}
	return false;
}

bool imgui_colormap_interface(rendering::shader_function::ColorMap::Uniforms& cm, const std::string& label)
{
	bool need_update = false;
	ImGui::TextUnformatted("Colormap:");
	ImGui::BeginGroup();
	need_update |= ImGui::RadioButton((std::string("BWR##") + label).c_str(), &(cm.color_map_), 0);
	ImGui::SameLine();
	need_update |= ImGui::RadioButton((std::string("CWR##") + label).c_str(), &cm.color_map_, 1);
	ImGui::SameLine();
	need_update |= ImGui::RadioButton((std::string("BCGYR##") + label).c_str(), &cm.color_map_, 2);
	ImGui::SameLine();
	need_update |= ImGui::RadioButton((std::string("BCR##") + label).c_str(), &cm.color_map_, 3);
	need_update |= ImGui::SliderInt((std::string("expansion##") + label).c_str(), &cm.expansion_, -2, 2);
	need_update |= ImGui::InputFloat((std::string("min##") + label).c_str(), &cm.min_value_);
	need_update |= ImGui::InputFloat((std::string("max##") + label).c_str(), &cm.max_value_);
	ImGui::EndGroup();
	return need_update;
}

} // namespace ui

} // namespace cgogn

#endif // CGOGN_UI_IMGUI_HELPERS_H_
