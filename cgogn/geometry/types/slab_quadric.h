#ifndef CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_
#define CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{
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
		Vec4 n1, n2;
		slab_normal(p1, p2, p3, n1, n2);
		_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
		_b = -2 * _A * (p1);
		_c = p1.transpose() * _A * p1;
	}


	bool slab_normal(const Vec4& p1, const Vec4& p2, const Vec4& p3, Vec4& n1, Vec4& n2)
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

		double dr12 = r2 - r1;
		double dr13 = r3 - r1;
		double dr23 = r3 - r2;

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
				contact_point2 = (r2 * c3 - r3 * c2) / (r2 - r3);
			}
			else if (dr23 < 1e-8)
			{
				contact_point1 = (r1 * c2 - r2 * c1) / (r1 - r2);
				contact_point2 = (r1 * c3 - r3 * c1) / (r1 - r3);
			}
			else
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}
			
			Vec3 c1cp1 = contact_point1 - c1;
			Vec3 c1cp2 = contact_point2 - c1;
			Vec3 cp1cp2 = contact_point2 - contact_point1;
			Vec3 inter_point;
			double dist;
			DistanceToLine(c1, contact_point1, contact_point2, dist, inter_point);
			double sangle = r1 / dist;
			if (fabs(sangle) > 1.)
				return false;
			double cangle = sqrt(1 - sangle * sangle);
			Vec3 c1inter = (inter_point - c1).normalized();
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
	
	Slab_Quadric operator+(const Slab_Quadric& q) const
	{
		Slab_Quadric sq;
		sq._A = _A + q._A;
		sq._b = _b + q._b;
		sq._c = _c + q._c;
		return sq;
	}

	inline Scalar eval(const Vec4& center) const
	{
		Scalar cost = center.transpose() * _A * center;
		cost += _b.transpose() * center;
		cost += _c;
		return cost;
		// return center.transpose() * _A * center + _b.transpose() * center + _c;
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

	void reset(){
		_A.setZero();
		_b.setZero();
		_c = 0;
	}

private:
	Mat4 _A;
	Vec4 _b;
	Scalar _c=0;
};
} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_