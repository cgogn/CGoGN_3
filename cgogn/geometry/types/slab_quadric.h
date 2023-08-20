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

 	inline Slab_Quadric(const Vec4& p1, const Vec4& p2, const Vec4& p3, const Vec4& n1, const Vec4&n2 )
	{
		Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		double r1 = p1[3];
		double r2 = p2[3];
		double r3 = p3[3];
		double d1 = (c2 - c1).norm();
		double d2 = (c3 - c1).norm();
		_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
		_b = -2 * _A * (p1);
		_c = p1.transpose() * _A * p1;
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
		cost -= _b.transpose() * center;
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

	void clear(){
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