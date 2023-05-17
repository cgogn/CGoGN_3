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

#ifndef CGOGN_GEOMETRY_TYPES_QUADRIC_H_
#define CGOGN_GEOMETRY_TYPES_QUADRIC_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

class Quadric
{
public:
	inline Quadric()
	{
		matrix_.setZero();
	}

	inline Quadric(const Vec3& p1, const Vec3& p2, const Vec3& p3)
	{
		Vec3 u = p2 - p1;
		Vec3 v = p3 - p1;
		Vec3 n = u.cross(v);
		n.normalize();
		Scalar d = -(p1.dot(n));
		Vec4 p = Vec4(n[0], n[1], n[2], d);
		matrix_ = p * p.transpose();
	}

	Quadric(const Quadric& q)
	{
		matrix_ = q.matrix_;
	}

	inline void zero()
	{
		matrix_.setZero();
	}

	Quadric& operator=(const Quadric& q)
	{
		matrix_ = q.matrix_;
		return *this;
	}

	Quadric& operator+=(const Quadric& q)
	{
		matrix_ += q.matrix_;
		return *this;
	}

	Scalar eval(const Vec3& v)
	{
		return eval(Vec4{v[0], v[1], v[2], 1.});
	}

	inline Scalar eval(const Vec4& v)
	{
		return v.transpose() * matrix_ * v;
	}

	bool optimized(Vec3& v)
	{
		Vec4 hv;
		bool b = optimized(hv);
		if (b)
		{
			v[0] = hv[0];
			v[1] = hv[1];
			v[2] = hv[2];
		}
		return b;
	}

	bool optimized(Vec4& v)
	{
		Mat4 m(matrix_);
		for (uint32 i = 0; i < 3; ++i)
			m(3, i) = 0.;
		m(3, 3) = 1.;
		Mat4 inverse;
		Scalar determinant;
		bool invertible;
		m.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			v = inverse * Vec4(0., 0., 0., 1.);
		return invertible;
	}

private:
	Mat4 matrix_;
};

class QMat_Quadric
{
public:
	inline Slab_Quadric()
	{
		A.setZero();
		b.setZero();
		c = 0;
	}
	
	

	inline Slab_Quadric(const Vec4& p1, const Vec4& p2, const Vec4& p3)
	{
		Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		double r1 = p1[3];
		double r2 = p2[3];
		double r3 = p3[3];
		double d1 = (c2 - c1).norm();
		double d2 = (c3 - c1).norm();
		auto [n1, n2] = Slab_normal(p1, p2, p3);
		A = n1.dot(n1) + n2.dot(n2);
		b = -2 * A * (p1);
	}

	inline Slab_Quadric(const Vec4& p1, const Vec4& p2)
	{
		Vec3 u = p2 - p1;
		Vec3 v = p3 - p1;
		Vec3 n = u.cross(v);
		n.normalize();
		Scalar d = -(p1.dot(n));
		Vec4 p = Vec4(n[0], n[1], n[2], d);
		matrix_ = p * p.transpose();
	}

	inline std::pair<Vec4,Vec4> Slab_normal(const Vec4& p1, const Vec4& p2, const Vec4& p3)
	{
		Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		double r1 = p1[3];
		double r2 = p2[3];
		double r3 = p3[3];
		
		double d1 = (c2 - c1).norm();
		double d2 = (c3 - c1).norm();
		
		Vec3 n1 = (c2 - c1).normalized();
		Vec3 n2 = (c3 - c1).normalized();
		Vec3 n = n1.cross(n2);

		double cosAlpha = (r2 - r1) / d1;
		double cosBeta = (r3 - r1) / d2;
		Vec3 p1 = c1 + (cosAlpha * r1) *n1;
		Vec3 p2 = c2 + (cosBeta * r2) * n2;

		Vec3 q1 = n.cross(n1);
		Vec3 q2 = n.cross(n2);

		double t = ((p2.x() - p1.x()) * q2.y() - (p2.y() - p1.y()) * q2.x() / (q1.x() * q2.y() - q1.y() * q2.x()));
		Vec3 p = p1 + t * q1;

		double rp = sqrt(r1*r1 - ((p - c1).norm() * (p-c1).norm()));
		Vec3 p_project_sphere = p + rp * n;
		Vec3 n_tangent  = (p_project_sphere - c1).normalized();
		Vec3 n_tangent_2 = ((p_project_sphere - c1) + (2 * (p1-p_project_sphere))).normalized();
		return {Vec4(n_tangent, 1), Vec4(n_tangent_2, 1)};
	}

	Quadric(const Quadric& q)
	{
		matrix_ = q.matrix_;
	}

	inline void zero()
	{
		matrix_.setZero();
	}

	Quadric& operator=(const Quadric& q)
	{
		matrix_ = q.matrix_;
		return *this;
	}

	Quadric& operator+=(const Quadric& q)
	{
		matrix_ += q.matrix_;
		return *this;
	}

	Scalar eval(const Vec3& v)
	{
		return eval(Vec4{v[0], v[1], v[2], 1.});
	}
		inline Scalar eval(const Vec4& v)
	{
		return v.transpose() * matrix_ * v;
	}

	bool optimized(Vec3& v)
	{
		Vec4 hv;
		bool b = optimized(hv);
		if (b)
		{
			v[0] = hv[0];
			v[1] = hv[1];
			v[2] = hv[2];
		}
		return b;
	}

	bool optimized(Vec4& v)
	{
		Mat4 m(matrix_);
		for (uint32 i = 0; i < 3; ++i)
			m(3, i) = 0.;
		m(3, 3) = 1.;
		Mat4 inverse;
		Scalar determinant;
		bool invertible;
		m.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			v = inverse * Vec4(0., 0., 0., 1.);
		return invertible;
	}

private:
	Mat4 A;
	Vec4 b;
	Scalar c;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_QUADRIC_H_
