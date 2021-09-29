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

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_QUADRIC_H_
