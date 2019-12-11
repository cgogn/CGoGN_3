#include <cgogn/core/functions/phi.h>

namespace cgogn {
Dart phi2bis(const MRCmap3& m,Dart d)
{
	uint32 face_id = m.cph().face_id(d);
	Dart it = d;

	it = phi2(m.mesh(),it);

	if(m.cph().face_id(it) == face_id)
		return it;
	else
	{
		do{
				it = phi2(m.mesh(),phi3(m.mesh(),it));
		}while(m.cph().face_id(it) != face_id);

		return it;
	}
}

Dart phi1(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.cph().dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	if(m.current_level() == m.cph().maximum_level())
		return phi1(m.mesh(),d);

	bool finished = false ;
	uint32 edge_id = m.cph().edge_id(d) ;
	Dart it = d ;
	do
	{
		it = phi1(m.mesh(),it) ;
		if(m.cph().dart_level(it) <= m.current_level())
			finished = true ;
		else
			while(m.cph().edge_id(it) != edge_id)
				it = phi1(m.mesh(),phi2bis(m,it)) ;
	} while(!finished) ;
	return it ;
}

Dart phi_1(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.cph().dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	bool finished = false ;
	Dart it = phi_1(m.mesh(),d) ;
	uint32 edge_id = m.cph().edge_id(it) ;

	do
	{
		if(m.cph().dart_level(it) <= m.current_level())
			finished = true ;
		else
		{
			it = phi_1(m.mesh(),it) ;
			while(m.cph().edge_id(it) != edge_id)
				it = phi_1(m.mesh(),phi2bis(m,it)) ;
		}
	} while(!finished) ;
	return it ;
}

Dart phi2(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.cph().dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	return phi2(m.mesh(),phi_1(m.mesh(),phi1(m,d)));
}

Dart phi3(const MRCmap3& m,Dart d)
{
	cgogn_message_assert(m.cph().dart_level(d) <= m.current_level(),
						 "Access to a dart introduced after current level") ;

	if(phi3(m.mesh(),d) == d)
		return d;

	return phi3(m.mesh(),phi_1(m.mesh(),phi1(m,d)));
}

}
