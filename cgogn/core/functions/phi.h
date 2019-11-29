#ifndef PHI_H
#define PHI_H

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn{

template<typename MAP>
Dart phi1(const MAP& m,Dart d){
	return (*(m.phi1_))[d.index];
}

template<typename MAP>
Dart phi_1(const MAP& m,Dart d)
{
	return (*(m.phi_1_))[d.index];
}

template<typename MAP>
Dart phi2(const MAP& m,Dart d)
{
	return (*(m.phi2_))[d.index];
}

template<typename MAP>
Dart phi3(const MAP& m,Dart d)
{
	return (*(m.phi3_))[d.index];
}

template <uint64 N,typename MAP>
inline Dart phi(const MAP& m,Dart d)
{
	static_assert(N % 10 <= mesh_traits<MAP>::dimension, "Composition of PHI: invalid index (phi1/phi2/phi3 only)");
	
	if constexpr(N%10 == 1 && mesh_traits<MAP>::dimension >= 1){
		return phi1(m,phi<N / 10>(m,d));
	}else{
		if constexpr(N%10 == 2 && mesh_traits<MAP>::dimension >= 2){
			return phi2(m,phi<N / 10>(m,d));
		}else{
			if constexpr(N%10 == 3 && mesh_traits<MAP>::dimension >= 3){
				return phi3(m,phi<N / 10>(m,d));
			}else{
				return d;
			}
		}
	}
}

template<typename MAP>
inline void phi1_sew(MAP& m, Dart d, Dart e)
{
	Dart f = phi1(m,d);
	Dart g = phi1(m,e);
	
	(*(m.phi1_))[d.index] = g;
	(*(m.phi1_))[e.index] = f;
	(*(m.phi_1_))[g.index] = d;
	(*(m.phi_1_))[f.index] = e;
}

template <typename MAP>
inline void phi1_unsew(MAP& m, Dart d)
{
	Dart e = phi1(m,d);
	Dart f = phi1(m,e);
	(*(m.phi1_))[d.index] = f;
	(*(m.phi1_))[e.index] = e;
	(*(m.phi_1_))[f.index] = d;
	(*(m.phi_1_))[e.index] = e;
}

template <typename MAP>
inline void phi2_sew(MAP& m,Dart d, Dart e)
{
	cgogn_assert(phi2(m,d) == d);
	cgogn_assert(phi2(m,e) == e);
	(*(m.phi2_))[d.index] = e;
	(*(m.phi2_))[e.index] = d;
}

template <typename MAP>
inline void phi2_unsew(MAP& m,Dart d)
{
	Dart e = phi2(m,d);
	(*(m.phi2_))[d.index] = d;
	(*(m.phi2_))[e.index] = e;
}

template <typename MAP>
inline void phi3_sew(MAP& m,Dart d, Dart e)
{
	cgogn_assert(phi3(m,d) == d);
	cgogn_assert(phi3(m,e) == e);
	(*(m.phi3_))[d.index] = e;
	(*(m.phi3_))[e.index] = d;
}

template <typename MAP>
inline void phi3_unsew(MAP& m,Dart d)
{
	Dart e = phi3(m,d);
	(*(m.phi3_))[d.index] = d;
	(*(m.phi3_))[e.index] = e;
}

}

#endif // PHI_H
