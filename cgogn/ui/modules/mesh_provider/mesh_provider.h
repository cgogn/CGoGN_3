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

#include <cgogn/ui/modules/mesh_provider/cgogn_module_mesh_provider_export.h>
#include <cgogn/ui/modules/mesh_provider/mesh_data.h>

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/vbo_update.h>

#include <string>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

class App;

class CGOGN_MODULE_MESH_PROVIDER_EXPORT MeshProvider : public Module
{
    using Mesh = CMap2;

public:

	MeshProvider(const App& app);
	~MeshProvider();

    Mesh* import_surface_from_file(const std::string& filename);

    template <typename FUNC>
    void foreach_mesh(const FUNC& f)
    {
        for (auto& [name, m] : meshes_)
            f(m.get(), name);
    }

    MeshData* mesh_data(const Mesh* m);

protected:

    void interface() override;

private:

    std::unordered_map<std::string, std::unique_ptr<Mesh>> meshes_;
    std::unordered_map<const Mesh*, MeshData> mesh_data_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MESH_PROVIDER_H_
