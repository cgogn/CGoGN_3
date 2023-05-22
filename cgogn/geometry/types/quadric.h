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

// Calculate the quadric matrix using the given points
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

 class Slab_Quadric
{
public:
	inline Slab_Quadric()
	{
		_A.setZero();
		_b.setZero();
		_c = 0;
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
 		auto [n1, n2] = slab_normal(p1, p2, p3);
		_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
 		_b = -2 * _A * (p1);
 		_c = p1.transpose() * _A * p1;
	}

	inline std::pair<Vec4, Vec4> slab_normal(const Vec4& p1, const Vec4& p2, const Vec4& p3)
	{
		 Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		double r1 = p1[3];
		double r2 = p2[3];
		double r3 = p3[3];
		// distance between two sphere
		double d1 = (c2 - c1).norm();
		double d2 = (c3 - c1).norm();
		// direction vector
		Vec3 n1 = (c2 - c1).normalized();
		Vec3 n2 = (c3 - c1).normalized();
		// normal vector of the triangle plane
		Vec3 n = n1.cross(n2);

		 double cosAlpha = (r2 - r1) / d1;
		double cosBeta = (r3 - r1) / d2;
		// projected point from sphere to the triangle plane
		Vec3 inter1 = c1 + (cosAlpha * r1) * n1;
		Vec3 inter2 = c2 + (cosBeta * r2) * n2;

		Vec3 q1 = n.cross(n1);
		Vec3 q2 = n.cross(n2);

		double t = ((inter2.x() - inter1.x()) * q2.y() - (inter2.y() - inter1.y()) * q2.x()) /
				   (q1.x() * q2.y() - q1.y() * q2.x());
		Vec3 p = inter1 + t * q1;

		double rp = sqrt(r1 * r1 - ((p - c1).norm() * (p - c1).norm()));
		Vec3 p_project_sphere = p + rp * n;
		Vec3 n_tangent = (p_project_sphere - c1).normalized();
		// reflect the tangent by the triangle plane
		Vec3 n_tangent_2 = ((p_project_sphere - c1) + (2 * (p - p_project_sphere))).normalized();
		 return {Vec4(n_tangent[0], n_tangent[1], n_tangent[2], 1), 
			Vec4(n_tangent_2[0], n_tangent_2[1], n_tangent_2[2], 1)};
	}
	 inline Slab_Quadric(const Vec4& p1, const Vec4& p2)
	{
		//Add virtual sphere to compute the tangent plan normal 
		Vec4 p3 = Vec4(0, 0, 0, (p1[3] + p2[3]) / 2);
		auto [n1, n2] = slab_normal(p1, p2, p3);
		Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		Eigen::AngleAxisd aa(0.5 * M_PI, (c2 - c1).normalized());
		Eigen::Quaterniond q(aa);
		Eigen::Quaterniond q_n1(n1), q_n2(n2);
		q_n1.w() = 0;
		q_n2.w() = 0;
		
		Eigen::Quaterniond rotated_n1 = q * q_n1 * q.inverse();
		Eigen::Quaterniond rotated_n2 = q * q_n2 * q.inverse();
		Vec4 r_n1 = Vec4(rotated_n1.x(),rotated_n1.y(),rotated_n1.z(),1);
		Vec4 r_n2 = Vec4(rotated_n2.x(), rotated_n2.y(), rotated_n2.z(), 1);
		_A = n1 * n1.transpose() + n2 * n2.transpose();
		_b = -2 * _A * (p1);
		_c = p1.transpose() * _A * p1;
		_A += r_n1 * r_n1.transpose() + r_n2 * r_n2.transpose();
		_b += -2 * _A * (p1);
		_c += p1.transpose() * _A * p1;
	}
	
	 Slab_Quadric& operator=(const Slab_Quadric& q)
	{
		_A = q._A;
		_b = q._b;
		_c = q._c;
		return *this;
	}

	Slab_Quadric& operator+=(const Slab_Quadric& q)
	{
		_A += q._A;
		_b += q._b;
		_c += q._c;
		return *this;
	}

	
	inline Scalar eval(const Vec4& center) const
	{
		Scalar v = center.transpose() * _A * center;
		v += _b.transpose() * center;
		v += _c;
		return v;
		//return center.transpose() * _A * center + _b.transpose() * center + _c;
	}

	bool optimized(Vec4& v)
	{
		Mat4 m(_A);
		Mat4 inverse;
		Scalar determinant;
		bool invertible;
		m.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			v = inverse * _b * 0.5;
		return invertible;
	}

private:
	Mat4 _A;
	Vec4 _b;
	Scalar _c;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_QUADRIC_H_
