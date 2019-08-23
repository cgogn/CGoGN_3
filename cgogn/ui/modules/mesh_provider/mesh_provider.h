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

#include <cgogn/io/surface_import.h>

#include <boost/synapse/emit.hpp>

#include <string>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

class App;

template <typename MESH>
class MeshProvider : public Module
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
    using AttributeGen = typename mesh_traits<MESH>::AttributeGen;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

public:

	MeshProvider(const App& app) :
		Module(app, "MeshProvider (" + mesh_traits<MESH>::name + ")"),
		show_mesh_inspector_(false),
		selected_mesh_(nullptr)
	{}
	~MeshProvider()
	{}

	MESH* load_surface_from_file(const std::string& filename)
	{
		std::string name = filename_from_path(filename);
		const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
		MESH* m = it->second.get();
		bool imported = cgogn::io::import_OFF(*m, filename);
		if (imported)
		{
			mesh_data_.emplace(m, MeshData(m));
			boost::synapse::emit<mesh_added>(this, m);
			return m;
		}
		else
		{
			meshes_.erase(name);
			return nullptr;
		}
	}

	template <typename FUNC>
	void foreach_mesh(const FUNC& f)
	{
		for (auto& [name, m] : meshes_)
			f(m.get(), name);
	}

	std::string mesh_name(const MESH* m)
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

	/////////////
	// SIGNALS //
	/////////////

	using mesh_added = struct mesh_added_ (*) (MESH* m);
	template <typename T>
	using attribute_changed_t = struct attribute_changed_t_(*)(Attribute<T>* attribute);
	using attribute_changed = struct attribute_changed_(*)(AttributeGen* attribute);
	using connectivity_changed = struct connectivity_changed_ (*) ();

	template <typename T>
	void emit_attribute_changed(const MESH* m, Attribute<T>* attribute)
	{
		MeshData<MESH>* md = mesh_data(m);
		md->update_vbo(attribute);
		if (static_cast<AttributeGen*>(md->bb_attribute_.get()) == static_cast<AttributeGen*>(attribute))
			md->update_bb();

		boost::synapse::emit<attribute_changed>(m, attribute);
		boost::synapse::emit<attribute_changed_t<T>>(m, attribute);
	}

	void emit_connectivity_changed(const MESH* m)
	{
		MeshData<MESH>* md = mesh_data(m);
		md->update_nb_cells();

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
				load_surface_from_file(result[0]);
			open_file_dialog = nullptr;
		}

		if (ImGui::BeginMenu(name_.c_str()))
        {
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, (bool)open_file_dialog);
			if (ImGui::MenuItem("Load mesh"))
				open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_files);
			ImGui::MenuItem("Show inspector", "", &show_mesh_inspector_);
			ImGui::PopItemFlag();
            ImGui::EndMenu();
        }
	}

	void interface() override
	{
		if (show_mesh_inspector_)
		{
			std::string name = mesh_traits<MESH>::name + " inspector";
			ImGui::Begin(name.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
			ImGui::SetWindowSize({0, 0});

			if (ImGui::ListBoxHeader("Select mesh"))
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
				ImGui::Text("Cells");
				ImGui::Columns(2);
				ImGui::Separator();
				ImGui::Text("Type"); ImGui::NextColumn();
				ImGui::Text("Number"); ImGui::NextColumn();
				ImGui::Separator();

				MeshData<MESH>* md = mesh_data(selected_mesh_);
				for (uint32 i = 0; i < std::tuple_size<typename mesh_traits<MESH>::Cells>::value; ++i)
				{
					ImGui::Text(mesh_traits<MESH>::CellNames[i].c_str()); ImGui::NextColumn();
					ImGui::Text("%d", md->nb_cells_[i]); ImGui::NextColumn();
				}
			}
		}
	}

private:

	std::vector<std::string> supported_files = { "Surface meshes", "*.off" };
	bool show_mesh_inspector_;
	const MESH* selected_mesh_;

	std::unordered_map<std::string, std::unique_ptr<MESH>> meshes_;
	std::unordered_map<const MESH*, MeshData<MESH>> mesh_data_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MESH_PROVIDER_H_
