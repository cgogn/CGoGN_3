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
		_add_A.setZero();
		_add_b.setZero();
		_add_c = 0;
	}

 	inline Slab_Quadric(const Vec4& p1, const Vec4& p2, const Vec4& p3, const Vec4& n1, const Vec4& n2,
						bool boundary = false)
	{
		if (!boundary)
		{
			_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
			_b = -2 * _A * p1;
			_c = p1.transpose() * _A * p1;
			_add_A.setZero();
			_add_b.setZero();
			_add_c = 0;
		}
		else
		{
			_A.setZero();
			_b.setZero();
			_c = 0;
			_add_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
			_add_b = -2 * _A * p1;
			_add_c = p1.transpose() * _A * p1;
		}
	}

	inline Slab_Quadric(const Vec4& p1, const Vec4& n1, double stability_ratio)
	{
		_A.setZero();
		_b.setZero();
		_c = 0;
		_add_A = n1 * n1.transpose() * 0.1 * stability_ratio*stability_ratio;
		_add_b = -2 * _A * p1;
		_add_c = p1.transpose() * _A * p1;
		
	}
	Slab_Quadric& operator=(const Slab_Quadric& q)
	{
		_A = q._A;
		_b = q._b;
		_c = q._c;
		_add_A = q._add_A;
		_add_b = q._add_b;
		_add_c = q._add_c;
		return *this;
	}

	Slab_Quadric& operator+=(const Slab_Quadric& q)
	{
		_A += q._A;
		_b += q._b;
		_c += q._c;
		_add_A += q._add_A;
		_add_b += q._add_b;
		_add_c += q._add_c;
		return *this;
	}
	
	friend Slab_Quadric operator+(Slab_Quadric lhs, Slab_Quadric& rhs)
	{
		lhs += rhs;
		return lhs;
	}

	inline Scalar eval(const Vec4& center) const
	{
		Scalar cost = center.transpose() * (_A + _add_A) * center;
		cost += (_b + _add_b).transpose() * center;
		cost += _c + _add_c;
		return cost;
		/*return center.transpose() * _A * center + _b.transpose() * center + _c;*/
	}

	bool optimized(Vec4& v)
	{
		Mat4 m(_A);
		Mat4 inverse;
		Scalar determinant;
		bool invertible;
		m.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			v = inverse * _b *0.5;
		return invertible;
	}

	void clear(){
		_A.setZero();
		_b.setZero();
		_c = 0;
		_add_A.setZero();
		_add_b.setZero();
		_add_c = 0;
 	}

private:
	Mat4 _A;
	Vec4 _b;
	Scalar _c;
	Mat4 _add_A;
	Vec4 _add_b;
	Scalar _add_c;
};
} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_