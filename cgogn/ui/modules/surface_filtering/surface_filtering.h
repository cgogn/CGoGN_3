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

#ifndef CGOGN_MODULE_SURFACE_FILTERING_H_
#define CGOGN_MODULE_SURFACE_FILTERING_H_

#include <cgogn/ui/modules/surface_filtering/cgogn_module_surface_filtering_export.h>

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace ui
{

class App;
class CMapProvider;

class CGOGN_MODULE_SURFACE_FILTERING_EXPORT SurfaceFiltering : public Module
{
    using Mesh = CMap2;

    template <typename T>
    using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;

    using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;

    using Vec3 = cgogn::geometry::Vec3;
    using Scalar = cgogn::geometry::Scalar;

public:

	SurfaceFiltering(const App& app);
	~SurfaceFiltering();

	void init();

	void filter_mesh(Mesh& m, Attribute<Vec3>* vertex_position);

protected:

    void interface() override;

private:

	Mesh* selected_mesh_;
	Attribute<Vec3>* selected_vertex_position_;
	CMapProvider* cmap_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_FILTERING_H_
