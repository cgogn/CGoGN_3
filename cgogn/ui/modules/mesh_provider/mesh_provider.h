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

#ifndef CGOGN_MODULE_MESH_PROVIDER_H_
#define CGOGN_MODULE_MESH_PROVIDER_H_

#include <cgogn/ui/modules/mesh_provider/mesh_data.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/portable-file-dialogs.h>

#include <cgogn/core/utils/string.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/io/graph/cg.h>
#include <cgogn/io/graph/skel.h>
#include <cgogn/io/surface/off.h>
#include <cgogn/io/volume/tet.h>

#include <boost/synapse/emit.hpp>

#include <string>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

class App;

template <typename MESH>
class MeshProvider : public ProviderModule
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
    using AttributeGen = typename mesh_traits<MESH>::AttributeGen;

    using Vertex = typename mesh_traits<MESH>::Vertex;

	using Scalar = geometry::Scalar;
    using Vec3 = geometry::Vec3;

public:

	MeshProvider(const App& app) :
		ProviderModule(app, "MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"),
		show_mesh_inspector_(false),
		selected_mesh_(nullptr),
		bb_min_(0, 0, 0),
		bb_max_(0, 0, 0)
	{}
	
	~MeshProvider()
	{}

	MESH* add_mesh(const std::string& name)
	{
		const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
		MESH* m = it->second.get();
		MeshData<MESH>& md = mesh_data_[m];
		md.init(m);
		boost::synapse::emit<mesh_added>(this, m);
		return m;
	}

	MESH* load_graph_from_file(const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 1)
		{
			std::string name = filename_from_path(filename);
			const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
			MESH* m = it->second.get();
			
			std::string ext = extension(filename);
			bool imported;
			if(ext.compare("cg") == 0)
				imported = cgogn::io::import_CG(*m, filename);
			else if(ext.compare("skel") == 0)
				imported = cgogn::io::import_SKEL(*m, filename);
			else
				imported = false;

			if (imported)
			{
				MeshData<MESH>& md = mesh_data_[m];
				md.init(m);
				std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
				if (vertex_position)
					set_mesh_bb_vertex_position(m, vertex_position);
				boost::synapse::emit<mesh_added>(this, m);
				return m;
			}
			else
			{
				meshes_.erase(name);
				return nullptr;
			}
		}
		else
			return nullptr;
	}

	MESH* load_surface_from_file(const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 2)
		{
			std::string name = filename_from_path(filename);
			const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
			MESH* m = it->second.get();
			bool imported = cgogn::io::import_OFF(*m, filename);
			if (imported)
			{
				MeshData<MESH>& md = mesh_data_[m];
				md.init(m);
				std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
				if (vertex_position)
					set_mesh_bb_vertex_position(m, vertex_position);
				boost::synapse::emit<mesh_added>(this, m);
				return m;
			}
			else
			{
				meshes_.erase(name);
				return nullptr;
			}
		}
		else
			return nullptr;
	}

	MESH* load_volume_from_file(const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 3)
		{
			std::string name = filename_from_path(filename);
			const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
			MESH* m = it->second.get();
			bool imported = cgogn::io::import_TET(*m, filename);
			if (imported)
			{
				MeshData<MESH>& md = mesh_data_[m];
				md.init(m);
				std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
				if (vertex_position)
					set_mesh_bb_vertex_position(m, vertex_position);
				boost::synapse::emit<mesh_added>(this, m);
				return m;
			}
			else
			{
				meshes_.erase(name);
				return nullptr;
			}
		}
		else
			return nullptr;
	}

	template <typename FUNC>
	void foreach_mesh(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH*>::value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter type");
		for (auto& [name, m] : meshes_)
			f(m.get(), name);
	}

	std::string mesh_name(const MESH* m) const
	{
		auto it = std::find_if(
			meshes_.begin(), meshes_.end(),
			[&] (const auto& pair) { return pair.second.get() == m; }
		);
		if (it != meshes_.end())
			return it->first;
		else
			return "";
	}

	MeshData<MESH>* mesh_data(const MESH* m)
	{
		auto it = mesh_data_.find(m);
		if (it != mesh_data_.end())
			return &(it->second);
		else
			return nullptr;
	}

	std::pair<Vec3, Vec3> meshes_bb() const override
	{
		return std::make_pair(bb_min_, bb_max_);
	}

private:

	void update_meshes_bb()
	{
		for (uint32 i = 0; i < 3; ++i)
		{
			bb_min_[i] = std::numeric_limits<float64>::max();
			bb_max_[i] = std::numeric_limits<float64>::lowest();
		}
		for (auto& [m, md] : mesh_data_)
		{
			for (uint32 i = 0; i < 3; ++i)
			{
				if (md.bb_min_[i] < bb_min_[i])
					bb_min_[i] = md.bb_min_[i];
				if (md.bb_max_[i] > bb_max_[i])
					bb_max_[i] = md.bb_max_[i];
			}
		}
	}

public:

	void set_mesh_bb_vertex_position(const MESH* m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		MeshData<MESH>& md = mesh_data_[m];
		md.bb_vertex_position_ = vertex_position;
		md.update_bb();
		update_meshes_bb();
		for (View* v : linked_views_)
			v->update_scene_bb();
	}

	/////////////
	// SIGNALS //
	/////////////

	using mesh_added = struct mesh_added_(*)(MESH* m);
	template <typename T>
	using attribute_changed_t = struct attribute_changed_t_(*)(Attribute<T>* attribute);
	using attribute_changed = struct attribute_changed_(*)(AttributeGen* attribute);
	using connectivity_changed = struct connectivity_changed_(*)();

	template <typename T>
	void emit_attribute_changed(const MESH* m, Attribute<T>* attribute)
	{
		MeshData<MESH>& md = mesh_data_[m];
		md.update_vbo(attribute);
		if (static_cast<AttributeGen*>(md.bb_vertex_position_.get()) == static_cast<AttributeGen*>(attribute))
		{
			md.update_bb();
			update_meshes_bb();
			for (View* v : linked_views_)
				v->update_scene_bb();
		}

		boost::synapse::emit<attribute_changed>(m, attribute);
		boost::synapse::emit<attribute_changed_t<T>>(m, attribute);
	}

	void emit_connectivity_changed(const MESH* m)
	{
		MeshData<MESH>* md = mesh_data(m);
		md->update_nb_cells();
		md->rebuild_cells_sets();
		md->set_primitives_dirty(rendering::POINTS);
		md->set_primitives_dirty(rendering::LINES);
		md->set_primitives_dirty(rendering::TRIANGLES);

		boost::synapse::emit<connectivity_changed>(m);
	}

protected:

	void main_menu() override
	{
		static std::shared_ptr<pfd::open_file> open_file_dialog;
		if (open_file_dialog && open_file_dialog->ready())
		{
			auto result = open_file_dialog->result();
			if (result.size())
			{
				if constexpr (mesh_traits<MESH>::dimension == 1)
					load_graph_from_file(result[0]);
				if constexpr (mesh_traits<MESH>::dimension == 2)
					load_surface_from_file(result[0]);
				if constexpr (mesh_traits<MESH>::dimension == 3)
					load_volume_from_file(result[0]);
			}
			open_file_dialog = nullptr;
		}

		if (ImGui::BeginMenu(name_.c_str()))
        {
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, (bool)open_file_dialog);
			if (ImGui::MenuItem("Load mesh"))
			{
				if constexpr (mesh_traits<MESH>::dimension == 1)
					open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_graph_files);
				if constexpr (mesh_traits<MESH>::dimension == 2)
					open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_surface_files);
				if constexpr (mesh_traits<MESH>::dimension == 3)
					open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_volume_files);
			}
			ImGui::MenuItem("Show inspector", "", &show_mesh_inspector_);
			ImGui::PopItemFlag();
            ImGui::EndMenu();
        }
	}

	void interface() override
	{
		if (show_mesh_inspector_)
		{
			std::string name = std::string{mesh_traits<MESH>::name} + " inspector";
			ImGui::Begin(name.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
			ImGui::SetWindowSize({0, 0});

			if (ImGui::ListBoxHeader("Mesh"))
			{
				foreach_mesh([this] (MESH* m, const std::string& name)
				{
					if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
						selected_mesh_ = m;
				});
				ImGui::ListBoxFooter();
			}

			if (selected_mesh_)
			{
				MeshData<MESH>* md = mesh_data(selected_mesh_);

				if (ImGui::BeginCombo("Position", md->bb_vertex_position_ ? md->bb_vertex_position_->name().c_str() : "-- select --"))
				{
					foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
					{
						bool is_selected = attribute == md->bb_vertex_position_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
							set_mesh_bb_vertex_position(selected_mesh_, attribute);
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
					ImGui::EndCombo();
				}
				
				ImGui::Separator();
				ImGui::TextUnformatted("Cells");
				ImGui::Columns(2);
				ImGui::Separator();
				ImGui::TextUnformatted("Type"); ImGui::NextColumn();
				ImGui::TextUnformatted("Number"); ImGui::NextColumn();
				ImGui::Separator();
				for (uint32 i = 0; i < std::tuple_size<typename mesh_traits<MESH>::Cells>::value; ++i)
				{
					ImGui::TextUnformatted(mesh_traits<MESH>::cell_names[i]); ImGui::NextColumn();
					ImGui::Text("%d", md->nb_cells_[i]); ImGui::NextColumn();
				}
			}
		}
	}

private:

	std::vector<std::string> supported_graph_files = { "Graph", "*.cg *.skel" };
	std::vector<std::string> supported_surface_files = { "Surface", "*.off" };
	std::vector<std::string> supported_volume_files = { "Volume", "*.tet" };

	bool show_mesh_inspector_;
	const MESH* selected_mesh_;

	std::unordered_map<std::string, std::unique_ptr<MESH>> meshes_;
	std::unordered_map<const MESH*, MeshData<MESH>> mesh_data_;
	Vec3 bb_min_, bb_max_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MESH_PROVIDER_H_
