/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_VOLUME_MR_MODELING_H_
#define CGOGN_MODULE_VOLUME_MR_MODELING_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/decimation/decimation.h>
#include <cgogn/modeling/algos/subdivision.h>

namespace cgogn
{

namespace ui
{

class VolumeMRModeling : public Module
{
	template <typename T>
	using Attribute = typename mesh_traits<CPH3>::template Attribute<T>;

	using Vertex = typename mesh_traits<CPH3>::Vertex;
	using Edge = typename mesh_traits<CPH3>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	VolumeMRModeling(const App& app)
		: Module(app, "VolumeMRModeling (" + std::string{mesh_traits<CPH3>::name} + ")"), selected_cph3_(nullptr),
		  selected_cmap3_(nullptr), selected_vertex_position_(nullptr)
	{
	}
	~VolumeMRModeling()
	{
	}

	CPH3* create_cph3(CPH3::CMAP& m, const std::string& name)
	{
		CPH3* cph = new CPH3(m);
		std::string cph_name;
		uint32 count = 1;
		do
		{
			cph_name = name + "_" + std::to_string(count);
			++count;
		} while (cph3_provider_->has_mesh(cph_name));
		cph3_provider_->register_mesh(cph, cph_name);
		return cph;
	}

	void subdivide(CPH3& m, Attribute<Vec3>* vertex_position)
	{
		uint32 cur = m.current_level_;
		m.current_level_ = m.maximum_level_;

		modeling::butterflySubdivisionVolumeAdaptative(m, 0.5f, vertex_position);

		m.current_level_ = cur;

		cph3_provider_->emit_connectivity_changed(&m);
		cph3_provider_->emit_attribute_changed(m, vertex_position);

		cmap3_provider_->emit_connectivity_changed(&static_cast<CPH3::CMAP&>(m));
		cmap3_provider_->emit_attribute_changed(&static_cast<CPH3::CMAP&>(m), vertex_position);
	}

protected:
	void init() override
	{
		cph3_provider_ = static_cast<ui::MeshProvider<CPH3>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<CPH3>::name} + ")"));

		cmap3_provider_ = static_cast<ui::MeshProvider<CPH3::CMAP>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<CPH3::CMAP>::name} + ")"));
	}

	void ui_interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("CMap3"))
		{
			cmap3_provider_->foreach_mesh([this](CPH3::CMAP* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_cmap3_))
				{
					selected_cmap3_ = m;
					selected_cmap3_name_ = name;
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_cmap3_)
		{
			if (ImGui::Button("Create CPH3"))
				create_cph3(*selected_cmap3_, selected_cmap3_name_);
		}

		if (ImGui::ListBoxHeader("CPH3"))
		{
			cph3_provider_->foreach_mesh([this](CPH3* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_cph3_))
				{
					selected_cph3_ = m;
					selected_vertex_position_.reset();
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_cph3_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			uint32 min = 0;
			if (ImGui::SliderScalar("Level", ImGuiDataType_U32, &selected_cph3_->current_level_, &min,
									&selected_cph3_->maximum_level_))
			{
				cph3_provider_->emit_connectivity_changed(selected_cph3_);
				cph3_provider_->emit_attribute_changed(selected_cph3_, selected_vertex_position_.get());
			}

			std::string selected_vertex_position_name_ =
				selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_cph3_,
												[this](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == selected_vertex_position_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														selected_vertex_position_ = attribute;
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (selected_vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##attribute"))
					selected_vertex_position_.reset();
			}

			if (selected_vertex_position_)
			{
				if (ImGui::Button("Subdivide"))
					subdivide(*selected_cph3_, selected_vertex_position_.get());
			}
		}

		ImGui::End();
	}

private:
	CPH3* selected_cph3_;
	CPH3::CMAP* selected_cmap3_;
	std::string selected_cmap3_name_;

	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;

	MeshProvider<CPH3>* cph3_provider_;
	MeshProvider<CPH3::CMAP>* cmap3_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_VOLUME_MR_MODELING_H_
