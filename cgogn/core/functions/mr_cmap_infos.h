#ifndef MR_CMAP_INFOS_H
#define MR_CMAP_INFOS_H

#include <cgogn/core/types/cmap/mr_cmap3.h>

namespace cgogn {

/**
* \brief add and init a new dart
* The dart is added to the current level of resolution
*/
Dart add_dart(MRCmap3& m);

inline Dart begin(const MRCmap3& m)
{ 
	Dart d = begin(m.mesh());
	while(d != end(m.mesh()) && m.cph().dart_level(d) > m.current_level())
		d = next(m.mesh(),d);
	return d; 
}
inline Dart next(const MRCmap3& m,Dart d)
{ 
	d = next(m.mesh(),d);
	while(d != end(m.mesh()) && m.cph().dart_level(d) > m.current_level())
		d = next(m.mesh(),d);
	return d;
}

template <typename CELL>
inline uint32 index_of(const MRCmap3& m,CELL c){
	static const Orbit orbit = CELL::ORBIT;
	Dart d = c.dart;
	
	if constexpr(orbit == PHI1_PHI2){
		d = volume_youngest_dart(m,c.dart);
	}
	if constexpr(orbit == PHI1_PHI3){
		d = face_youngest_dart(m,c.dart);
	}
	if constexpr(orbit == PHI2_PHI3){
		d = edge_youngest_dart(m,c.dart);
	}
	
	return index_of(m.mesh(),CELL(d));
}

/**
 * Return the level of the edge of d in the current level map
 * \details The level of an edge is the maximum of the levels of
 * its darts. As phi1(d) and phi2(d) are from the same level we can
 * optimize by checking phi1(d) instead of phi2(d)
 */
uint32 edge_level(const MRCmap3& m, Dart d);

/**
 * Return the youngest dart of the edge of d in the current level map
 */
Dart edge_youngest_dart(const MRCmap3& m, Dart d);

/**
 * Return the level of the face of d in the current level map
 * \details The level of a face is the minimum of the levels of its edges
 * but a specific treatment has to be done in the particular case of a
 * face with all neighboring faces are regularly subdivided
 * but not the face itself
 */
uint32 face_level(const MRCmap3 &m, Dart d);

/**
 * Given the face of d in the current level map,
 * return a level 0 dart of its origin face
 */
Dart face_origin(const MRCmap3 &m, Dart d);

/**
 * Return the oldest dart of the face of d in the current level map
 */
Dart face_oldest_dart(const MRCmap3& m, Dart d);

/**
 * Return the youngest dart of the face of d in the current level map
 */
Dart face_youngest_dart(const MRCmap3& m, Dart d);


/**
 * Return the level of the volume of d in the current level map
 * \details The level of a volume is the minimum of the levels of its faces
 * but a specific treatment has to be done in the particular case of a
 * volume with all neighboring volumes are regularly subdivided
 * but not the volume itself
 */
uint32 volume_level(const MRCmap3 &m, Dart d);

/**
 * Return the oldest dart of the volume of d in the current level map
 */
Dart volume_oldest_dart(const MRCmap3 &m, Dart d);

/**
 * Return the youngest dart of the volume of d in the current level map
 */
Dart volume_youngest_dart(const MRCmap3 &m, Dart d);

/**
 * Return true if the edge of d in the current level map
 * has already been subdivided to the next level
 * As before testing phi2(d) or phi1(d) is the same
 */
bool edge_is_subdivided(const MRCmap3& m, Dart d);

/**
 * Return true if the face of d in the current level map
 * has already been subdivided to the next level
 */
bool face_is_subdivided(const MRCmap3 &m, Dart d);

/**
 * Return true if the face of d in the current level map
 * is subdivided to the next level
 * and none of its resulting faces is in turn subdivided to the next level
 * \details
 * A face whose level in the current level map is lower than the current
 * level can not be subdivided to higher levels
 */
bool face_is_subdivided_once(const MRCmap3 &m, Dart d);

/**
 * Return true if the volume of d in the current level map
 * has already been subdivided to the next level
 */
bool volume_is_subdivided(const MRCmap3 &m, Dart d);

}

#endif // MR_CMAP_INFOS_H
