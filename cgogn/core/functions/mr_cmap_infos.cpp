#include <cgogn/core/functions/mr_cmap_infos.h>
#include <cgogn/core/functions/phi.h>
#include <cgogn/core/functions/traversals/volume.h>


namespace cgogn{

uint32 edge_level(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");
	return std::max(m.dart_level(d),m.dart_level(phi1(m,d)));
}

Dart edge_youngest_dart(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");
	if(m.dart_level(d) > m.dart_level(phi2(m,d)))
		return d;
	return phi2(m,d);
}

uint32 face_level(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	if(m.current_level() == 0)
		return 0;

	Dart it = d;
	Dart old = it;
	uint32 l_old = m.dart_level(old);
	uint32 fLevel = edge_level(m,it);
	do
	{
		it = phi1(m,it);
		uint32 dl = m.dart_level(it);

		// compute the oldest dart of the face in the same time
		if(dl < l_old)
		{
			old = it;
			l_old = dl;
		}
		uint32 l = edge_level(m,it);
		fLevel = l < fLevel ? l : fLevel;
	} while(it != d);
	MRCmap3 m2(m,fLevel);

	uint32 nbSubd = 0;
	it = old;
	uint32 eId = m2.edge_id(old);
	unsigned int init_dart_level = m2.dart_level(it);
	do
	{
		++nbSubd;
		it = phi1(m2,it);
	} while(m2.edge_id(it) == eId && m2.dart_level(it) != init_dart_level);

	while(nbSubd > 1)
	{
		nbSubd /= 2;
		--fLevel;
	}
	return fLevel ;
}

Dart face_origin(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");
	Dart p = d;
	uint32 pLevel = m.dart_level(p);
	MRCmap3 m2(m,pLevel);
	while(pLevel > 0){
		p = face_oldest_dart(m2,p);
		pLevel = m2.dart_level(p);
		m2 = std::move(MRCmap3(m,pLevel));
	}
	return p;
}

Dart face_oldest_dart(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	Dart it = d ;
	Dart oldest = it ;
	uint32 l_old = m.dart_level(oldest);
	do
	{
		uint32 l = m.dart_level(it);
		if(l == 0)
			return it;

		if(l < l_old)
		//		if(l < l_old || (l == l_old && it < oldest))
		{
			oldest = it;
			l_old = l;
		}
		it = phi1(m,it);
	} while(it != d);

	return oldest;
}

Dart face_youngest_dart(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	Dart it = d ;
	Dart youngest = it ;
	uint32 l_young = m.dart_level(youngest);
	do
	{
		uint32 l = m.dart_level(it);
		if(l == m.current_level())
			return it;

		if(l > l_young)
		//		if(l < l_young || (l == l_young && it < youngest))
		{
			youngest = it;
			l_young = l;
		}
		it = phi1(m,it);
	} while(it != d);

	return youngest;
}

uint32 volume_level(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	//The level of a volume is the
	// minimum of the levels of its faces

	Dart oldest = d ;
	uint32 lold = m.dart_level(oldest);
	uint32 vLevel = std::numeric_limits<uint32>::max() ;

	foreach_incident_face(m,MRCmap3::Volume(d),[&](MRCmap3::Face d)->bool{
		uint32 fLevel = face_level(m,d.dart);
		vLevel = fLevel < vLevel ? fLevel : vLevel;
		Dart old = face_oldest_dart(m,d.dart);
		if(m.dart_level(old) < lold){
			oldest = old;
			lold = m.dart_level(old);
		}
		return true;
	});

	MRCmap3 m2(m,vLevel);

	uint32 nbSubd = 0;
	Dart it = oldest;
	uint32 eId = m2.edge_id(oldest);
	do{
		++nbSubd;
		it = phi<121>(m2,it);
	}while(m2.edge_id(it) == eId && lold != m2.dart_level(it));

	while(nbSubd > 1){
		nbSubd /= 2;
		--vLevel;
	}

	return vLevel;
}

Dart volume_oldest_dart(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	Dart oldest = d;
	uint32 l_old = m.dart_level(oldest);
	foreach_incident_face(m,MRCmap3::Volume(oldest), [&](MRCmap3::Face f)->bool{
		Dart old = face_oldest_dart(m,f.dart);
		uint32 l = m.dart_level(old);
		if(l < l_old){
			oldest = old;
			l_old = l;
		}
		return true;
	});

	return oldest;
}

Dart volume_youngest_dart(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	Dart youngest = d;
	uint32 l_young = m.dart_level(youngest);
	foreach_incident_face(m,MRCmap3::Volume(youngest), [&](MRCmap3::Face f)->bool{
		Dart young = face_youngest_dart(m,f.dart);
		uint32 l = m.dart_level(young);
		if(l > l_young)
		{
			youngest = young;
			l_young = l;
		}
		return true;
	});

	return youngest;
}

/**
 * Return true if the edge of d in the current level map
 * has already been subdivided to the next level
 * As before testing phi2(d) or phi1(d) is the same
 */
bool edge_is_subdivided(const MRCmap3 &m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");
	
	if(m.current_level() == m.maximum_level())
		return false ;
	
	Dart d1 = phi1(m,d);
	MRCmap3 m2(m,m.current_level()+1);
	Dart d1_l = phi1(m2,d);
	if(d1 != d1_l)
		return true;
	else
		return false;
}

/**
 * Return true if the face of d in the current level map
 * has already been subdivided to the next level
 */
bool face_is_subdivided(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");
	
	uint32 fLevel = face_level(m,d) ;
	if(fLevel < m.current_level())
		return false ;
	
	bool subd = false ;
	
	MRCmap3 m2(m,m.current_level()+1);
	if(m2.dart_level(phi1(m2,d)) == m2.current_level()
		&& m2.edge_id(phi1(m2,d)) != m2.edge_id(d))
		subd = true;
	
	return subd;
}

/**
 * Return true if the face of d in the current level map
 * is subdivided to the next level
 * and none of its resulting faces is in turn subdivided to the next level
 * \details
 * A face whose level in the current level map is lower than the current
 * level can not be subdivided to higher levels
 */
bool face_is_subdivided_once(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");

	uint32 fLevel = face_level(m,d);
	if(fLevel < m.current_level())
		return false;

	uint32 degree = 0 ;
	bool subd = false ;
	bool subdOnce = true ;
	Dart fit = d ;
	MRCmap3 m2(m,m.current_level()+1),m3(m,m.current_level()+2);
	do
	{
		if(m2.dart_level(phi1(m2,fit)) == m2.current_level()
			&& m2.edge_id(phi1(m2,fit)) != m2.edge_id(fit))
		{
			subd = true ;
			if(m3.dart_level(phi1(m3,fit)) == m3.current_level()
				&& m3.edge_id(phi1(m3,fit)) != m3.edge_id(fit))
				subdOnce = false ;
		}
		++degree ;
		fit = phi1(m,fit) ;
	} while(subd && subdOnce && fit != d) ;

	if(degree == 3 && subd)
	{
		Dart cf = phi2(m2,phi1(m2,d)) ;
		if(m3.dart_level(phi1(m3,cf)) == m3.current_level()
			&& m3.edge_id(phi1(m3,cf)) != m3.edge_id(cf))
		subdOnce = false ;
	}

	return subd && subdOnce ;
}

bool volume_is_subdivided(const MRCmap3& m, Dart d){
	cgogn_message_assert(m.dart_level(d) <= m.current_level(),
				 "Access to a dart introduced after current level");
	
	uint vLevel = volume_level(m,d);
	if(vLevel < m.current_level())
		return false ;
	bool faceAreSubdivided = face_is_subdivided(m,d);
	
	
	foreach_incident_face(m,MRCmap3::Volume(d),[&](MRCmap3::Face f)->bool{
		faceAreSubdivided &= face_is_subdivided(m,f.dart);
		return true;
	});
	
	bool subd = false;
	MRCmap3 m2(m,m.current_level()+1);
	if(faceAreSubdivided && m2.dart_level(phi<112>(m2,d)) == m2.current_level() && m2.face_id(phi<112>(m2,d)) != m2.face_id(d))
		subd = true;
		
	return subd;
}


}
