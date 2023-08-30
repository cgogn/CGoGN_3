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

#ifndef CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_
#define CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>

#include <cgogn/geometry/types/slab_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

template <typename NONMANIFOLD>
struct DecimationSQEM_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	using Vertex = typename NONMANIFOLD::Vertex;
	using Edge = typename NONMANIFOLD::Edge;
	using Face = typename NONMANIFOLD::Face;

	using Vec3 = typename geometry::Vec3;
	using Vec4 = typename geometry::Vec4;
	using Scalar = typename geometry::Scalar;
	using Slab_Quadric = typename geometry::Slab_Quadric;

	DecimationSQEM_Helper(NONMANIFOLD& m, std::shared_ptr<Attribute<Vec4>>& sphere_info,
						  std::shared_ptr<Attribute<Slab_Quadric>>& sphere_quadric,
		std::shared_ptr<Attribute<std::pair<Vec4,Vec4>>>& slab_normals,
		std::shared_ptr<Attribute<double>>& stability_ratio)
		: m_(m), sphere_info_(sphere_info), sphere_quadric_(sphere_quadric), slab_normals_(slab_normals),
		  stability_ratio_(stability_ratio)
	{
		
	}

	~DecimationSQEM_Helper()
	{
	}

	void initial_slab_mesh()
	{
		parallel_foreach_cell(m_, [&](Vertex v) -> bool {
			value<Slab_Quadric>(m_, sphere_quadric_, v).clear();
			return true;
		});
		foreach_cell(m_, [&](Face f) -> bool {
			Vec4 n1 = Vec4(0, 0, 0, 1), n2 = Vec4(0, 0, 0, 1);
			auto ifv = incident_vertices(m_, f);
			Vec4 p1 = value<Vec4>(m_, sphere_info_, ifv[0]);
			Vec4 p2 = value<Vec4>(m_, sphere_info_, ifv[1]);
			Vec4 p3 = value<Vec4>(m_, sphere_info_, ifv[2]);
			if (slab_normal(p1, p2, p3, n1, n2) && n1 != Vec4(0, 0, 0, 1) && n2 != Vec4(0, 0, 0, 1))
			{
				value<std::pair<Vec4, Vec4>>(m_, slab_normals_, f) = {n1, n2};
				Slab_Quadric q1(p1, p2, p3, n1, n2);
				Slab_Quadric q2(p2, p1, p3, n1, n2);
				Slab_Quadric q3(p3, p1, p2, n1, n2);
				value<Slab_Quadric>(m_, sphere_quadric_, ifv[0]) += q1;
				value<Slab_Quadric>(m_, sphere_quadric_, ifv[1]) += q2;
				value<Slab_Quadric>(m_, sphere_quadric_, ifv[2]) += q3;
			}
			return true;
		});
	}

	void initial_boundary_mesh()
	{
		foreach_cell(m_, [&](Edge e) {
			if (is_boundary(m_, e))
			{
				auto ief = incident_faces(m_, e);
				auto iev = incident_vertices(m_, e);
				auto iv1_e = incident_edges(m_, iev[0]);
				auto iv2_e = incident_edges(m_, iev[1]);
				double sr = value<double>(m_, stability_ratio_, e);
				Vec4 p1 = value<Vec4>(m_, sphere_info_, iev[0]);
				Vec3 p2 = value<Vec4>(m_, sphere_info_, iev[1]);
				Vec3 v1 = p1.head<3>();
				Vec3 v2 = p2.head<3>();
				Vec3 v1v2 = v2 - v1;
				Vec3 r1 = p1.w();
				Vec3 r2 = p2.w();
				if (ief.size() == 1)
				{
					auto iefv = incident_vertices(m_, ief);

					for (auto vert : iefv)
					{
						auto vert_pos = value<Vec4>(m_, sphere_info_, vert).head<3>();
						if (vert_pos != v1 && vert_pos != v2)
						{
							v3 = vert_pos;
							break;
						}
					}
					Vec3 v1v3 = v3 - v1;
					Vec3 iefv_nor = v1v2.cross(v1v3);
					Vec3 normal_boundary = iefv_nor.cross(v1v2);
					double temp_cos = acos(normal_boundary.dot(v1v3) / v1v3.norm() / normal_boundary.norm());
					bool dir = temp_cos > M_PI / 2.0 ? true : false;
					if (!dir)
						normal_boundary *= -1;
					Vec4 n1 = Vec4(normal_boundary.x(), normal_boundary.y(), normal_boundary.z(), 1.0);
					Slab_Quadric boundary_v1(v1, n1, sr);
					Slab_Quadric boundary_v2(v2, n1, sr);
					value<Slab_Quadric>(m_, sphere_quadric_, iev[0]) += boundary_v1;
					value<Slab_Quadric>(m_, sphere_quadric_, iev[1]) += boundary_v2;
				}
				else if (ief.size() == 0)
				{
					p3 = Vec4(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z(), (r1 + r2) * 0.5);
					Vec4 n1 = Vec4(0, 0, 0, 1), n2 = Vec4(0, 0, 0, 1);
					if (slab_normal(p1, p2, p3, n1, n2))
					{
						if (n1 != Vec4(0, 0, 0, 1) && n2 != Vec4(0, 0, 0, 1))
						{
							// Add slab quadric of current slab
							Slab_Quadric q1(p1, p2, p3, n1, n2, true);
							Slab_Quadric q2(p2, p1, p3, n1, n2, true);
							value<Slab_Quadric>(m_, sphere_quadric_, iev[0]) += q1;
							value<Slab_Quadric>(m_, sphere_quadric_, iev[1]) += q2;
							Vec3 t1 = v1v2.cross(n1);
							Vec3 t2 = v1v2.cross(n2);
							// Add slab quadric of two orthogonal plane of the slab
							Vec4 tnormal1 = Vec4(t1.x(), t1.y(), t1.z(), 1.0);
							Vec4 tnormal2 = Vec4(t2.x(), t2.y(), t2.z(), 1.0);
							Slab_Quadric q3(p1, p2, p3, tnormal1, tnormal2, true);
							Slab_Quadric q4(p2, p1, p3, tnormal1, tnormal2, true);
							value<Slab_Quadric>(m_, sphere_quadric_, iev[0]) += q3;
							value<Slab_Quadric>(m_, sphere_quadric_, iev[1]) += q4;
						}
					}
					if (iv1_e.size() == 1)
					{
						Vec4 boundary_vertex_norm = Vec4(v1v2.x(), v1v2.y(), v1v2.z(), 1.0);
						Slab_Quadric boundary_v1(v1, boundary_vertex_norm, sr);
						value<Slab_Quadric>(m_, sphere_quadric_, iev[0]) += boundary_v1;
					}
					if (iv1_e.size() == 1)
					{
						Vec4 boundary_vertex_norm = Vec4(-v1v2.x(), -v1v2.y(), -v1v2.z(), 1.0);
						Slab_Quadric boundary_v2(v2, boundary_vertex_norm, sr);
						value<Slab_Quadric>(m_, sphere_quadric_, iev[1]) += boundary_v2;
					}
				}
			}
			return true;
			});
	}

	bool slab_normal(Vec4& p1, Vec4& p2, Vec4& p3, Vec4& n1, Vec4& n2)
	{

		Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		double r1 = p1[3];
		double r2 = p2[3];
		double r3 = p3[3];

		Vec3 c12 = c2 - c1;
		Vec3 c13 = c3 - c1;
		Vec3 c23 = c3 - c2;
		double dc12 = c12.norm();
		double dc13 = c13.norm();
		double dc23 = c23.norm();

		double dr12 = fabs(r2 - r1);
		double dr13 = fabs(r3 - r1);
		double dr23 = fabs(r3 - r2);

		if ((dc12 < 1e-8) || (dc13 < 1e-8) || (dc23 < 1e-8))
			return false;

		Vec3 n = c12.cross(c13).normalized();

		if ((dr12 < 1e-8) && (dr13 < 1e-8) && (dr23 < 1e-8))
		{
			n1 = Vec4(n[0], n[1], n[2], 1);
			n2 = Vec4(-n[0], -n[1], -n[2], 1);
			return true;
		}
		else
		{
			// two points on the tangent plane
			Vec3 contact_point1, contact_point2;

			if (dr12 < 1e-8)
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}
			else if (dr13 < 1e-8)
			{
				contact_point1 = (r2 * c1 - r1 * c2) / (r2 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}
			else if (dr23 < 1e-8)
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r1 * c2 - r2 * c1) / (r1 - r2);
			}
			else
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}

			Vec3 inter_point;
			double dist;
			DistanceToLine(c1, contact_point1, contact_point2, dist, inter_point);
			double sangle = r1 / dist;
			if (fabs(sangle) > 1.)
				return false;
			double cangle = sqrt(1. - sangle * sangle);
			Vec3 c1inter = (c1 - inter_point ).normalized();
			Vec3 n1_vec3 = n * cangle - c1inter * sangle;
			Vec3 n2_vec3 = -n * cangle - c1inter * sangle;
			n1_vec3.normalize();
			n2_vec3.normalize();
			Vec3 v11 = c1 + r1 * n1_vec3;
			Vec3 v12 = c2 + r2 * n1_vec3;
			Vec3 v13 = c3 + r3 * n1_vec3;

			Vec3 nv1 = (v12- v11).cross(v13-v11).normalized();

			Vec3 v21 = c1 + r1 * n2_vec3;
			Vec3 v22 = c2 + r2 * n2_vec3;
			Vec3 v23 = c3 + r3 * n2_vec3;

			Vec3 nv2 = (v22 - v21).cross(v23 - v21).normalized() * -1;

			n1 = Vec4(nv1[0], nv1[1], nv1[2], 1);
			n2 = Vec4(nv2[0], nv2[1], nv2[2], 1);
		}
		return true;
	}

	bool DistanceToLine(const Vec3& p, const Vec3& v0, const Vec3& v1, double& dist, Vec3& fp)
	{
		Vec3 v0v1(v1 - v0), pv0(v0 - p), pv1(v1 - p);
		double area = fabs(v0v1.cross(pv0).norm());
		if (v0v1.norm() > 1e-12)
		{
			dist = area / v0v1.norm();
			double t = (pv0.dot(pv0) - pv0.dot(pv1)) / (pv0.dot(pv0) + pv1.dot(pv1) - 2 * pv0.dot(pv1));
			fp = (1 - t) * v0 + t * v1;
			return true;
		}
		else
			return false;
	}

	Scalar edge_cost(Edge e, const Vec4& p)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		Slab_Quadric q;
		Slab_Quadric& q0 = value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
		Slab_Quadric& q1 = value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
		q+= value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
		q+= value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
		return q.eval(p);
	}

	Vec4 edge_optimal(Edge e)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		Slab_Quadric q;
		q += value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
		q += value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
		Vec4 p;
		if (q.optimized(p))
			if (p[3] < 0)
				return Scalar(0.5) * (value<Vec4>(m_, sphere_info_, iv[0]) + value<Vec4>(m_, sphere_info_, iv[1]));
			else
				return p;
		else
			return Scalar(0.5) * (value<Vec4>(m_, sphere_info_, iv[0]) + value<Vec4>(m_, sphere_info_, iv[1]));
	}


	NONMANIFOLD& m_;
	std::shared_ptr<Attribute<Vec4>> sphere_info_;
	std::shared_ptr<Attribute<Slab_Quadric>> sphere_quadric_;
	std::shared_ptr <Attribute<std::pair<Vec4, Vec4>>> slab_normals_;
	std::shared_ptr<Attribute<double>> stability_ratio_;

};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_
