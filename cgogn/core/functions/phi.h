#ifndef PHI_H
#define PHI_H

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn{

template<typename MAP>
Dart phi1(const MAP& m,Dart d){
	return (*(m.phi1_))[d.index];
}

inline Dart phi1(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	if(m.current_level() == m.maximum_level())
		return phi1(m.mesh(),d);

	bool finished = false ;
	uint32 edge_id = m.edge_id(d) ;
	Dart it = d ;
	do
	{
		it = phi1(m.mesh(),it) ;
		if(m.dart_level(it) <= m.current_level())
			finished = true ;
		else
			while(m.edge_id(it) != edge_id)
				it = phi1(m.mesh(),phi2bis(m,it)) ;
	} while(!finished) ;
	return it ;
}

template<typename MAP>
Dart phi_1(const MAP& m,Dart d)
{
	return (*(m.phi_1_))[d.index];
}

inline Dart phi_1(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	bool finished = false ;
	Dart it = phi_1(m.mesh(),d) ;
	uint32 edge_id = m.edge_id(it) ;

	do
	{
		if(m.dart_level(it) <= m.current_level())
			finished = true ;
		else
		{
			it = phi_1(m.mesh(),it) ;
			while(m.edge_id(it) != edge_id)
				it = phi_1(m.mesh(),phi2bis(m,it)) ;
		}
	} while(!finished) ;
	return it ;
}

template<typename MAP>
Dart phi2(const MAP& m,Dart d)
{
	return (*(m.phi2_))[d.index];
}

inline Dart phi2(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	return phi2(m.mesh(),phi_1(m.mesh(),phi1(m,d)));
}

template<typename MAP>
Dart phi3(const MAP& m,Dart d)
{
	return (*(m.phi3_))[d.index];
}

inline Dart phi3(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	if(phi3(m.mesh(),d) == d)
		return d;

	return phi3(m.mesh(),phi_1(m.mesh(),phi1(m,d)));
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

inline Dart phi2bis(const MRCmap3& m,Dart d)
{
	uint32 face_id = m.face_id(d);
	Dart it = d;

	it = phi2(m.mesh(),it);

	if(m.face_id(it) == face_id)
		return it;
	else
	{
		do{
				it = phi2(m.mesh(),phi3(m.mesh(),it));
		}while(m.face_id(it) != face_id);

		return it;
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
