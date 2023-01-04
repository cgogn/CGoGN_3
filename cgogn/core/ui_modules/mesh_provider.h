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

#ifndef CGOGN_MODULE_MESH_PROVIDER_H_
#define CGOGN_MODULE_MESH_PROVIDER_H_

#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/portable-file-dialogs.h>

#include <cgogn/core/functions/mesh_ops/global.h>
#include <cgogn/core/ui_modules/mesh_data.h>
#include <cgogn/core/utils/string.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/io/graph/cg.h>
#include <cgogn/io/graph/cgr.h>
#include <cgogn/io/graph/skel.h>
#include <cgogn/io/incidence_graph/ig.h>
#include <cgogn/io/surface/obj.h>
#include <cgogn/io/surface/off.h>
#include <cgogn/io/surface/ply.h>
// #include <cgogn/io/volume/cgns.h>
#include <cgogn/io/volume/mesh.h>
#include <cgogn/io/volume/meshb.h>
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
	MeshProvider(const App& app)
		: ProviderModule(app, "MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  bb_min_(0, 0, 0), bb_max_(0, 0, 0)
	{
		// for (auto& n : new_attribute_name_)
		// 	n[0] = '\0';
		if constexpr (mesh_traits<MESH>::dimension == 1)
			supported_formats_ = &supported_graph_formats_;
		if constexpr (mesh_traits<MESH>::dimension == 2)
			supported_formats_ = &supported_surface_formats_;
		if constexpr (mesh_traits<MESH>::dimension == 3)
			supported_formats_ = &supported_volume_formats_;
	}

	~MeshProvider()
	{
	}

	MESH* add_mesh(const std::string& name)
	{
		if constexpr (std::is_default_constructible_v<MESH>)
		{
			const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
			MESH* m = it->second.get();
			if (inserted)
			{
				MeshData<MESH>& md = mesh_data(*m);
				md.init(m);
				boost::synapse::emit<mesh_added>(this, m);
			}
			return m;
		}
		else
			return nullptr;
	}

	MESH* clone_mesh(const MESH& m)
	{
		const std::string& m_name = mesh_name(m);
		std::string name =
			remove_extension(m_name) + "_" + std::to_string(number_of_meshes()) + "." + extension(m_name);
		MESH* result = add_mesh(name);
		copy(*result, m);
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*result, "position");
		if (vertex_position)
			set_mesh_bb_vertex_position(*result, vertex_position);
		boost::synapse::emit<mesh_added>(this, result);
		emit_connectivity_changed(*result);
		// TODO: emit attributes changed ?
		return result;
	}

	void register_mesh(MESH* m, const std::string& name)
	{
		const auto [it, inserted] = meshes_.emplace(name, std::unique_ptr<MESH>(m));
		if (inserted)
		{
			MeshData<MESH>& md = mesh_data_[m];
			md.init(m);
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_mesh_bb_vertex_position(*m, vertex_position);
			boost::synapse::emit<mesh_added>(this, m);
		}
	}

	void clear_mesh(MESH& m)
	{
		clear(m);
		emit_connectivity_changed(m);
		// TODO: emit attributes changed ?
	}

	void copy_mesh(MESH& dst, const MESH& src)
	{
		copy(dst, src);
		emit_connectivity_changed(dst);
	}

	void remove_mesh(MESH& m)
	{
		// TODO
	}

	bool has_mesh(const std::string& name) const
	{
		return meshes_.count(name) == 1;
	}

	MESH* load_graph_from_file(const std::string& filename)
	{
		std::string name = filename_from_path(filename);
		if (has_mesh(name))
			name = remove_extension(name) + "_" + std::to_string(number_of_meshes()) + "." + extension(name);
		const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
		MESH* m = it->second.get();

		std::string ext = extension(filename);
		bool imported;

		if constexpr (mesh_traits<MESH>::dimension == 1 && std::is_default_constructible_v<MESH>)
		{
			if (ext.compare("cg") == 0)
				imported = cgogn::io::import_CG(*m, filename);
			else if (ext.compare("cgr") == 0)
				imported = cgogn::io::import_CGR(*m, filename);
			else if (ext.compare("ig") == 0)
				imported = cgogn::io::import_IG(*m, filename);
			else if (ext.compare("skel") == 0)
				imported = cgogn::io::import_SKEL(*m, filename);
			else
				imported = false;
		}
		else if constexpr (std::is_same_v<MESH, IncidenceGraph>)
		{
			if (ext.compare("cg") == 0)
				imported = cgogn::io::import_CG(*m, filename);
			else if (ext.compare("ig") == 0)
				imported = cgogn::io::import_IG(*m, filename);
			else
				imported = false;
		}

		if (imported)
		{
			MeshData<MESH>& md = mesh_data(*m);
			md.init(m);
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_mesh_bb_vertex_position(*m, vertex_position);
			boost::synapse::emit<mesh_added>(this, m);
			return m;
		}
		else
		{
			meshes_.erase(name);
			return nullptr;
		}
	}

	void save_graph_to_file(MESH& m, const Attribute<Vec3>* vertex_position, const std::string& filetype,
							const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 1)
		{
			if (filetype.compare("cg") == 0)
				cgogn::io::export_CG(m, vertex_position, filename + ".cg");
			else if (filetype.compare("ig") == 0)
				cgogn::io::export_IG(m, vertex_position, filename + ".ig");
			// else if (filetype.compare("cgr") == 0)
			// 	// TODO cgogn::io::export_CGR();
			// else if (filetype.compare("skel") == 0)
			// 	// TODO cgogn::io::export_SKEL();
		}
	}

	MESH* load_surface_from_file(const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 2 && std::is_default_constructible_v<MESH>)
		{
			std::string name = filename_from_path(filename);
			if (has_mesh(name))
				name = remove_extension(name) + "_" + std::to_string(number_of_meshes()) + "." + extension(name);
			const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
			MESH* m = it->second.get();

			std::string ext = extension(filename);
			bool imported = false;
			if (ext.compare("off") == 0)
				imported = cgogn::io::import_OFF(*m, filename);
			else if (ext.compare("obj") == 0)
				imported = cgogn::io::import_OBJ(*m, filename);
			else if (ext.compare("ply") == 0)
				imported = cgogn::io::import_PLY(*m, filename);
			else if (ext.compare("ig") == 0)
			{
				if constexpr (std::is_same_v<MESH, IncidenceGraph>)
					imported = cgogn::io::import_IG(*m, filename);
			}

			if (imported)
			{
				MeshData<MESH>& md = mesh_data(*m);
				md.init(m);
				std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
				if (vertex_position)
					set_mesh_bb_vertex_position(*m, vertex_position);
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

	void save_surface_to_file(MESH& m, const Attribute<Vec3>* vertex_position, const std::string& filetype,
							  const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 2)
		{
			if (filetype.compare("off") == 0)
				cgogn::io::export_OFF(m, vertex_position, filename + ".off");
			else if (filetype.compare("ig") == 0)
				cgogn::io::export_IG(m, vertex_position, filename + ".ig");
		}
	}

	MESH* load_volume_from_file(const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 3 && std::is_default_constructible_v<MESH>)
		{
			std::string name = filename_from_path(filename);
			if (has_mesh(name))
				name = remove_extension(name) + "_" + std::to_string(number_of_meshes()) + "." + extension(name);
			const auto [it, inserted] = meshes_.emplace(name, std::make_unique<MESH>());
			MESH* m = it->second.get();

			std::string ext = extension(filename);
			bool imported;
			if (ext.compare("tet") == 0)
				imported = cgogn::io::import_TET(*m, filename);
			else if (ext.compare("mesh") == 0 || ext.compare("meshb") == 0)
				imported = cgogn::io::import_MESHB(*m, filename);
			else
				imported = false;

			if (imported)
			{
				MeshData<MESH>& md = mesh_data(*m);
				md.init(m);
				std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
				if (vertex_position)
					set_mesh_bb_vertex_position(*m, vertex_position);
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

	void save_volume_to_file(MESH& m, const Attribute<Vec3>* vertex_position, const std::string& filetype,
							 const std::string& filename)
	{
		if constexpr (mesh_traits<MESH>::dimension == 3)
		{
			if (filetype.compare("mesh") == 0)
				cgogn::io::export_MESH(m, vertex_position, filename + ".mesh");
			// else if (filetype.compare("cgns") == 0)
			// 	cgogn::io::export_CGNS(m, vertex_position, filename + ".cgns");

			// else if (filetype.compare("tet") == 0)
			// 	// TODO cgogn::io::export_TET();
			// else if (filetype.compare("meshb") == 0)
			// 	// TODO cgogn::io::export_MESHB();
		}
	}

	template <typename FUNC>
	void foreach_mesh(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH&>::value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter type");
		for (auto& [name, m] : meshes_)
			f(*m, name);
	}

	inline uint32 number_of_meshes()
	{
		return uint32(meshes_.size());
	}

	std::string mesh_name(const MESH& m) const
	{
		auto it =
			std::find_if(meshes_.begin(), meshes_.end(), [&](const auto& pair) { return pair.second.get() == &m; });
		if (it != meshes_.end())
			return it->first;
		else
			return "";
	}

	MeshData<MESH>& mesh_data(const MESH& m)
	{
		return mesh_data_[&m];
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
	void set_mesh_bb_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		MeshData<MESH>& md = mesh_data(m);
		md.bb_vertex_position_ = vertex_position;
		md.update_bb();
		update_meshes_bb();
		for (View* v : linked_views_)
		{
			v->update_scene_bb();
			v->show_entire_scene();
		}
	}

	/////////////
	// SIGNALS //
	/////////////

	using mesh_added = struct mesh_added_ (*)(MESH* m);
	template <typename T>
	using attribute_changed_t = struct attribute_changed_t_ (*)(Attribute<T>* attribute);
	using attribute_changed = struct attribute_changed_ (*)(AttributeGen* attribute);
	using connectivity_changed = struct connectivity_changed_ (*)();
	template <typename CELL>
	using cells_set_changed = struct cells_set_changed_ (*)(CellsSet<MESH, CELL>* set);

	template <typename T>
	void emit_attribute_changed(const MESH& m, Attribute<T>* attribute)
	{
		MeshData<MESH>& md = mesh_data(m);
		md.update_vbo(attribute);
		if (static_cast<AttributeGen*>(md.bb_vertex_position_.get()) == static_cast<AttributeGen*>(attribute))
		{
			md.update_bb();
			update_meshes_bb();
			for (View* v : linked_views_)
				v->update_scene_bb();
		}

		for (View* v : linked_views_)
			v->request_update();

		boost::synapse::emit<attribute_changed>(&m, attribute);
		boost::synapse::emit<attribute_changed_t<T>>(&m, attribute);
	}

	void emit_connectivity_changed(const MESH& m)
	{
		MeshData<MESH>& md = mesh_data(m);
		md.update_nb_cells();
		md.rebuild_cells_sets();
		md.set_all_primitives_dirty();

		for (View* v : linked_views_)
			v->request_update();

		boost::synapse::emit<connectivity_changed>(&m);
	}

	template <typename CELL>
	void emit_cells_set_changed(const MESH& m, CellsSet<MESH, CELL>* set)
	{
		boost::synapse::emit<cells_set_changed<CELL>>(&m, set);
	}

protected:
	void main_menu() override
	{
		static std::shared_ptr<pfd::open_file> open_file_dialog;
		if (open_file_dialog && open_file_dialog->ready())
		{
			auto result = open_file_dialog->result();
			if (uint32(result.size()) > 0)
			{
				if constexpr (mesh_traits<MESH>::dimension == 1)
				{
					for (auto file : result)
						load_graph_from_file(file);
				}
				if constexpr (mesh_traits<MESH>::dimension == 2)
				{
					for (auto file : result)
						load_surface_from_file(file);
				}
				if constexpr (mesh_traits<MESH>::dimension == 3)
				{
					for (auto file : result)
						load_volume_from_file(file);
				}
			}
			open_file_dialog = nullptr;
		}

		open_save_popup_ = false;
		if (ImGui::BeginMenu(name_.c_str()))
		{
			if (ImGui::MenuItem("Add mesh"))
				add_mesh(std::string{mesh_traits<MESH>::name});
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, (bool)open_file_dialog);
			if (ImGui::MenuItem("Load mesh"))
			{
				if constexpr (mesh_traits<MESH>::dimension == 1)
					open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_graph_files_,
																		pfd::opt::multiselect);
				if constexpr (mesh_traits<MESH>::dimension == 2)
					open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_surface_files_,
																		pfd::opt::multiselect);
				if constexpr (mesh_traits<MESH>::dimension == 3)
					open_file_dialog = std::make_shared<pfd::open_file>("Choose file", ".", supported_volume_files_,
																		pfd::opt::multiselect);
			}
			ImGui::PopItemFlag();
			if (ImGui::MenuItem("Save mesh"))
				open_save_popup_ = true;
			ImGui::EndMenu();
		}
	}

	void popups() override
	{
		if (open_save_popup_)
			ImGui::OpenPopup("Save");

		if (ImGui::BeginPopupModal("Save", NULL, ImGuiWindowFlags_AlwaysAutoResize))
		{
			static MESH* selected_mesh = nullptr;
			static char filename[32] = "\0";
			static std::string filetype = (*supported_formats_)[0];
			static std::function<void()> cleanup = []() {};
			bool close_popup = false;

			imgui_mesh_selector(this, selected_mesh, "Mesh", [&](MESH& m) { selected_mesh = &m; });
			if (ImGui::BeginCombo("Filetype", filetype.c_str()))
			{
				for (const std::string& t : *supported_formats_)
				{
					bool is_selected = t == filetype;
					if (ImGui::Selectable(t.c_str(), is_selected))
						filetype = t;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
			ImGui::InputText("Filename", filename, 32);

			if (selected_mesh)
			{
				static std::shared_ptr<Attribute<Vec3>> selected_vertex_position = nullptr;
				imgui_combo_attribute<Vertex, Vec3>(
					*selected_mesh, selected_vertex_position, "Position",
					[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position = attribute; });
				if (selected_vertex_position)
				{
					if (ImGui::Button("Save", ImVec2(120, 0)))
					{
						if constexpr (mesh_traits<MESH>::dimension == 1)
							save_graph_to_file(*selected_mesh, selected_vertex_position.get(), filetype, filename);
						if constexpr (mesh_traits<MESH>::dimension == 2)
							save_surface_to_file(*selected_mesh, selected_vertex_position.get(), filetype, filename);
						if constexpr (mesh_traits<MESH>::dimension == 3)
							save_volume_to_file(*selected_mesh, selected_vertex_position.get(), filetype, filename);
						close_popup = true;
					}
				}
				cleanup = [&]() { selected_vertex_position = nullptr; };
			}

			if (ImGui::Button("Cancel", ImVec2(120, 0)))
				close_popup = true;

			if (close_popup)
			{
				selected_mesh = nullptr;
				filename[0] = '\0';
				cleanup();
				ImGui::CloseCurrentPopup();
			}
			ImGui::EndPopup();
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(this, selected_mesh_, "Mesh", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			MeshData<MESH>& md = mesh_data(*selected_mesh_);

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, md.bb_vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_mesh_bb_vertex_position(*selected_mesh_, attribute);
												});

			if (ImGui::Button("Clone"))
				clone_mesh(*selected_mesh_);

			ImGui::Separator();
			ImGui::TextUnformatted("Size");
			ImGui::Separator();

			if (ImGui::BeginTable("MeshSize", 2))
			{
				ImGui::TableSetupColumn("CellType");
				ImGui::TableSetupColumn("Number");
				ImGui::TableHeadersRow();

				for (uint32 i = 0; i < std::tuple_size<typename mesh_traits<MESH>::Cells>::value; ++i)
				{
					ImGui::TableNextColumn();
					ImGui::TextUnformatted(mesh_traits<MESH>::cell_names[i]);
					ImGui::TableNextColumn();
					ImGui::Text("%d", md.nb_cells_[i]);
				}
				ImGui::EndTable();
			}

			ImGui::Separator();
			ImGui::TextUnformatted("Attributes");
			ImGui::Separator();

			if (ImGui::BeginTable("MeshAttributes", 2))
			{
				ImGui::TableSetupColumn("CellType");
				ImGui::TableSetupColumn("Name");
				ImGui::TableHeadersRow();

				auto names = md.attributes_names();
				for (uint32 i = 0; i < std::tuple_size<typename mesh_traits<MESH>::Cells>::value; ++i)
				{
					ImGui::TableNextColumn();
					ImGui::TextUnformatted(mesh_traits<MESH>::cell_names[i]);
					ImGui::TableNextColumn();
					ImGui::PushItemWidth(-1);
					if (ImGui::ListBoxHeader((std::string("##") + mesh_traits<MESH>::cell_names[i]).c_str(),
											 names[i].size()))
					{
						for (auto& n : names[i])
							ImGui::Text("%s", n.c_str());
						ImGui::ListBoxFooter();
					}
					// ImGui::PopItemWidth();
					// ImGui::NextColumn();
					// ImGui::NextColumn();
					// ImGui::PushItemWidth(-1);
					// ImGui::InputText((std::string("##") + mesh_traits<MESH>::cell_names[i]).c_str(),
					// new_attribute_name_[i], 				 32); ImGui::PopItemWidth(); ImGui::SameLine(); if
					// (ImGui::Button((std::string("Add##") + mesh_traits<MESH>::cell_names[i]).c_str()))
					// {
					// }
				}
				ImGui::EndTable();
			}
		}
	}

private:
	std::vector<std::string> supported_graph_formats_ = {"cg", "ig", "cgr", "skel"};
	std::vector<std::string> supported_graph_files_ = {"Graph", "*.cg *.ig *.cgr *.skel"};

	std::vector<std::string> supported_surface_formats_ = {"off", "obj", "ply", "ig"};
	std::vector<std::string> supported_surface_files_ = {"Surface", "*.off *.obj *.ply *.ig"};

	std::vector<std::string> supported_volume_formats_ = {"mesh"};			 //, "cgns"};
	std::vector<std::string> supported_volume_files_ = {"Volume", "*.mesh"}; // *.cgns"};

	std::vector<std::string>* supported_formats_ = nullptr;

	bool open_save_popup_ = false;

	const MESH* selected_mesh_;
	// std::array<char[32], std::tuple_size<typename mesh_traits<MESH>::Cells>::value> new_attribute_name_;

	std::unordered_map<std::string, std::unique_ptr<MESH>> meshes_;
	std::unordered_map<const MESH*, MeshData<MESH>> mesh_data_;
	Vec3 bb_min_, bb_max_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MESH_PROVIDER_H_
