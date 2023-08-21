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
						  std::shared_ptr<Attribute<Slab_Quadric>> sphere_quadric
						  )
		: m_(m), sphere_info_(sphere_info), sphere_quadric_(sphere_quadric)
	{
		parallel_foreach_cell(m_, [&](Vertex v) -> bool {
			value<Slab_Quadric>(m_, sphere_quadric_, v).clear();
			return true;
		});
		foreach_cell(m_, [&](Face f) -> bool {
			Vec4 n1, n2;
			if (slab_normal(f, n1, n2))
			{
				std::vector<Vertex> iv = incident_vertices(m_, f);
				Slab_Quadric q1(value<Vec4>(m_, sphere_info_, iv[0]), value<Vec4>(m_, sphere_info_, iv[1]),
								value<Vec4>(m_, sphere_info_, iv[2]), n1, n2);
				Slab_Quadric q2(value<Vec4>(m_, sphere_info_, iv[1]), value<Vec4>(m_, sphere_info_, iv[0]),
								value<Vec4>(m_, sphere_info_, iv[2]), n1, n2);
				Slab_Quadric q3(value<Vec4>(m_, sphere_info_, iv[2]), value<Vec4>(m_, sphere_info_, iv[0]),
								value<Vec4>(m_, sphere_info_, iv[1]), n1, n2);
				value<Slab_Quadric>(m_, sphere_quadric_, iv[0]) += q1;
				value<Slab_Quadric>(m_, sphere_quadric_, iv[1]) += q2;
				value<Slab_Quadric>(m_, sphere_quadric_, iv[2]) += q3;
				
			}
			return true;
		});
	}

	~DecimationSQEM_Helper()
	{
	}

	bool slab_normal(Face f, Vec4& n1, Vec4& n2)
	{
		std::vector<Vertex> iv = incident_vertices(m_, f);
		
		Vec4 p1 = value<Vec4>(m_, sphere_info_, iv[0]);
		Vec4 p2 = value<Vec4>(m_, sphere_info_, iv[1]);
		Vec4 p3 = value<Vec4>(m_, sphere_info_, iv[2]);

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
			n1 = Vec4(n1_vec3[0], n1_vec3[1], n1_vec3[2], 1);
			n2 = Vec4(n2_vec3[0], n2_vec3[1], n2_vec3[2], 1);
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


};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_
