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

#ifndef CGOGN_GEOMETRY_TYPES_SPHEREICAL_QUADRIC_H_
#define CGOGN_GEOMETRY_TYPES_SPHEREICAL_QUADRIC_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

struct Spherical_Quadric
{
	Spherical_Quadric()
	{
		this->clear();
	}

	void clear()
	{
		_A.setZero();
		_b.setZero();
		_c = 0;
	}

	inline Spherical_Quadric(const Vec4& p, const Vec4& n)
	{
		_A = 2 * (n * n.transpose());
		_b = 2 * (n.transpose() * p * n.transpose());
		_c = n.transpose() * p * p.transpose() * n;
	}

	Spherical_Quadric& operator=(const Spherical_Quadric& q)
	{
		_A = q._A;
		_b = q._b;
		_c = q._c;
		return *this;
	}

	Spherical_Quadric& operator+=(const Spherical_Quadric& q)
	{
		_A += q._A;
		_b += q._b;
		_c += q._c;
		return *this;
	}

	Spherical_Quadric& operator*(Scalar s)
	{
		_A *= s;
		_b *= s;
		_c *= s;
		return *this;
	}

	friend Spherical_Quadric operator+(const Spherical_Quadric& lhs, const Spherical_Quadric& rhs)
	{
		Spherical_Quadric result(lhs);
		result += rhs;
		return result;
	}

	Scalar eval(const Vec4& p) const
	{
		return (0.5 * p.transpose() * _A * p) - (Scalar)(_b.transpose() * p) + _c;
	}

	bool optimized(Vec4& sphere)
	{
		Mat4 inverse;
		inverse.setZero();
		Scalar determinant;
		bool invertible;
		_A.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			sphere = inverse * _b;
		return invertible;
	}

	friend std::ostream& operator<<(std::ostream& os, const Spherical_Quadric& quad)
	{
		os << quad._A << ", " << quad._b << ", " << quad._c;
		return os;
	}

	Mat4 _A;
	Vec4 _b;
	Scalar _c = 0;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_SPHERICAL_QUADRIC_H_
