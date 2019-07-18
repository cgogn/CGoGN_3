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

#ifndef CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_H_
#define CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_H_

#include <cgogn/ui/modules/surface_differential_properties/cgogn_module_surface_differential_properties_export.h>

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace ui
{

class App;
class MeshProvider;

class CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_EXPORT SurfaceDifferentialProperties : public Module
{
    using Mesh = CMap2;

    template <typename T>
    using Attribute = typename mesh_traits<Mesh>::Attribute<T>;

    using Vertex = typename mesh_traits<Mesh>::Vertex;
    using Edge = typename mesh_traits<Mesh>::Edge;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	SurfaceDifferentialProperties(const App& app);
	~SurfaceDifferentialProperties();
    
	void init();

	void compute_normal(
		const Mesh& m,
		const Attribute<Vec3>* vertex_position,
		Attribute<Vec3>* vertex_normal
	);

	void compute_curvature(
		const Mesh& m,
		Scalar radius,
		const Attribute<Vec3>* vertex_position,
		const Attribute<Vec3>* vertex_normal,
		const Attribute<Scalar>* edge_angle,
		Attribute<Scalar>* vertex_kmax,
		Attribute<Scalar>* vertex_kmin,
		Attribute<Vec3>* vertex_Kmax,
		Attribute<Vec3>* vertex_Kmin,
		Attribute<Vec3>* vertex_Knormal
	);

protected:

    void interface() override;

private:

	Mesh* selected_mesh_;
	const Attribute<Vec3>* selected_vertex_position_;
	Attribute<Vec3>* selected_vertex_normal_;
	Attribute<Scalar>* selected_vertex_kmax_;
	Attribute<Scalar>* selected_vertex_kmin_;
	Attribute<Vec3>* selected_vertex_Kmax_;
	Attribute<Vec3>* selected_vertex_Kmin_;
	Attribute<Vec3>* selected_vertex_Knormal_;
	MeshProvider* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_H_
