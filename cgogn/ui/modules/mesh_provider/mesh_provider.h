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

public:

	MeshProvider(const App& app) :
		Module(app, "MeshProvider (" + mesh_traits<MESH>::name + ")")
	{}
	~MeshProvider()
	{}

	MESH* import_surface_from_file(const std::string& filename)
	{
		const auto [it, inserted] = meshes_.emplace(filename_from_path(filename), std::make_unique<MESH>());
		MESH* m = it->second.get();
		cgogn::io::import_OFF(*m, filename);
		mesh_data_.emplace(m, MeshData(m));

		boost::synapse::emit<mesh_added>(this, m);

		return m;
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
	using attribute_changed_t = struct attribute_changed_(*)(Attribute<T>* attribute);
	using attribute_changed = struct attribute_changed_(*)(AttributeGen* attribute);
	using connectivity_changed = struct connectivity_changed_ (*) ();

	template <typename T>
	void emit_attribute_changed(const MESH* m, Attribute<T>* attribute)
	{
		boost::synapse::emit<attribute_changed>(m, attribute);
		boost::synapse::emit<attribute_changed_t<T>>(m, attribute);
	}

protected:

	void interface() override
	{}

private:

	std::unordered_map<std::string, std::unique_ptr<MESH>> meshes_;
	std::unordered_map<const MESH*, MeshData<MESH>> mesh_data_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MESH_PROVIDER_H_
