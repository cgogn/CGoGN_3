/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
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

#ifndef CGOGN_GEOMETRY_DISTANCE_HEAT_SOLVER_H_
#define CGOGN_GEOMETRY_DISTANCE_HEAT_SOLVER_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/cells_set.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
class DistanceHeatSolver
{

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using Edge = typename mesh_traits<MESH>::Edge;

private:
	Eigen::SparseMatrix<Scalar> massMatrix;
	// Eigen::SparseMatrix Lc;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>> scalarHeatSolver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>> poissonSolver;
	double shortTime;
	uint32 nb_vertices;

	std::shared_ptr<Attribute<uint32>> vertex_index;
	std::shared_ptr<Attribute<Vec3>> position;
	std::shared_ptr<Attribute<Vec3>> face_normal;
	MESH* _mesh;

public:
	DistanceHeatSolver(MESH& mesh, std::shared_ptr<Attribute<Vec3>> position,
					   std::shared_ptr<Attribute<Vec3>> face_normal, double time_multiplier = 1.0)
	{
		this->_mesh = &mesh;
		this->position = position;
		this->face_normal = face_normal;
		vertex_index = get_or_add_attribute<uint32, Vertex>(mesh, "__distance_heat_index");

		nb_vertices = 0;
		foreach_cell(mesh, [&](Vertex v) -> bool {
			value<uint32>(mesh, vertex_index, v) = nb_vertices++;
			return true;
		});

		auto area = get_or_add_attribute<Scalar, Vertex>(mesh, "__area");
		geometry::compute_area<Vertex>(mesh, position.get(), area.get());

		auto Lc = geometry::cotan_operator_matrix(mesh, vertex_index.get(), position.get());

		Eigen::VectorXd areas = Eigen::VectorXd(nb_vertices);

		foreach_cell(mesh, [&](Vertex v) -> bool {
			areas(value<uint32>(mesh, vertex_index, v)) = value<Scalar>(mesh, area, v);
			return true;
		});

		double t = 0;
		int nb_edges = 0;
		foreach_cell(mesh, [&](Edge e) -> bool {
			t += geometry::length(mesh, e, position.get());
			nb_edges++;
			return true;
		});

		t /= nb_edges;
		t *= t;
		t *= time_multiplier;
		shortTime = t;

		Eigen::SparseMatrix<Scalar> heatOp =
			Eigen::SparseMatrix<Scalar, Eigen::ColMajor>(areas.asDiagonal()) - (t * Lc);
		scalarHeatSolver.compute(heatOp);
		poissonSolver.compute(Lc);
		remove_attribute<Vertex>(mesh, area);
	}

	~DistanceHeatSolver()
	{
		remove_attribute<Vertex>(*_mesh, vertex_index);
	}

	void solve(Attribute<Scalar>* out_heat_distance, ui::CellsSet<MESH, Vertex>& cell_source)
	{
		Eigen::VectorXd source = Eigen::VectorXd::Zero(nb_vertices);

		cell_source.foreach_cell([&](Vertex vec) { source(value<uint32>(*_mesh, vertex_index.get(), vec)) = 1; });

		auto u = scalarHeatSolver.solve(source);

		auto heat_attribute = get_or_add_attribute<Scalar, Vertex>(*_mesh, "__heat");
		parallel_foreach_cell(*_mesh, [&](Vertex v) -> bool {
			value<Scalar>(*_mesh, heat_attribute, v) = u(value<uint32>(*_mesh, vertex_index, v));
			return true;
		});

		// gradient
		// auto gradient_attribute = get_or_add_attribute<Vec3, Face>(*_mesh, "__heat_gradient");
		auto gradient_normalized_attribute = get_or_add_attribute<Vec3, Face>(*_mesh, "__heat_normalized_gradient");
		foreach_cell(*_mesh, [&](Face f) -> bool {
			Vec3 Ni = value<Vec3>(*_mesh, face_normal, f);
			std::vector<Vertex> vertices = incident_vertices(*_mesh, f);
			Vec3 grad;
			for (int i = 0; i < 3; i++)
			{
				Scalar Ui = value<Scalar>(*_mesh, heat_attribute, vertices[(i + 2) % 3]);
				Vec3 ei =
					value<Vec3>(*_mesh, position, vertices[(i + 1) % 3]) - value<Vec3>(*_mesh, position, vertices[i]);
				grad += Ni.cross(ei) * Ui;
			}
			grad /= 2.0 * geometry::area(*_mesh, f, position.get());
			// value<Vec3>(*_mesh, gradient_attribute, f) = grad;
			value<Vec3>(*_mesh, gradient_normalized_attribute, f) = (grad / grad.norm() * -1.0);
			return true;
		});

		auto integrated_divergence = get_or_add_attribute<Scalar, Vertex>(*_mesh, "__integrated_divergence");
		foreach_cell(*_mesh, [&](Vertex v) -> bool {
			Vec3 vpos = value<Vec3>(*_mesh, position, v);
			Scalar div = 0.0;
			foreach_incident_face(*_mesh, v, [&](Face f) -> bool {
				Vec3 Xj = value<Vec3>(*_mesh, gradient_normalized_attribute, f);
				std::vector<Vertex> vertices = incident_vertices(*_mesh, f);
				std::remove_if(vertices.begin(), vertices.end(), [&](Vertex item) -> bool {
					return value<uint32_t>(*_mesh, vertex_index, item) == value<uint32_t>(*_mesh, vertex_index, v);
				});
				Vec3 p0 = value<Vec3>(*_mesh, position, vertices[0]);
				Vec3 p1 = value<Vec3>(*_mesh, position, vertices[1]);
				Vec3 e1 = p0 - vpos;
				Vec3 e2 = p1 - vpos;
				double teta1 = geometry::angle(-e2, p0 - p1);
				double teta2 = geometry::angle(-e1, p1 - p0);
				div += 1.0 / tan(teta1) * e1.dot(Xj) + 1.0 / tan(teta2) * e2.dot(Xj);

				return true;
			});
			value<Scalar>(*_mesh, integrated_divergence, v) = div * 0.5;
			return true;
		});

		Eigen::VectorXd b = Eigen::VectorXd::Zero(nb_vertices);
		foreach_cell(*_mesh, [&](Vertex v) -> bool {
			b(value<uint32>(*_mesh, vertex_index.get(), v)) = value<Scalar>(*_mesh, integrated_divergence, v);
			return true;
		});

		auto phi = poissonSolver.solve(b);

		parallel_foreach_cell(*_mesh, [&](Vertex v) -> bool {
			value<Scalar>(*_mesh, out_heat_distance, v) = phi(value<uint32>(*_mesh, vertex_index.get(), v));
			return true;
		});

		Scalar min = std::numeric_limits<Scalar>::max();
		foreach_cell(*_mesh, [&](Vertex v) -> bool {
			if (value<Scalar>(*_mesh, out_heat_distance, v) < min)
			{
				min = value<Scalar>(*_mesh, out_heat_distance, v);
			}
			return true;
		});
		parallel_foreach_cell(*_mesh, [&](Vertex v) -> bool {
			value<Scalar>(*_mesh, out_heat_distance, v) -= min;
			return true;
		});

		// remove_attribute<Face>(*_mesh, gradient_attribute);
		remove_attribute<Face>(*_mesh, gradient_normalized_attribute);
		remove_attribute<Vertex>(*_mesh, integrated_divergence);
		remove_attribute<Vertex>(*_mesh, heat_attribute);
	}
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_DISTANCE_HEAT_SOLVER_H_
