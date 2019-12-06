#ifndef MR_CMAP_INFOS_H
#define MR_CMAP_INFOS_H

#include <cgogn/core/types/cmap/mr_cmap3.h>

namespace cgogn {

/**
* \brief add and init a new dart
* The dart is added to the current level of resolution
*/
inline Dart add_dart(MRCmap3& m)
{
	Dart d = add_dart(m.mesh());
	m.inc_nb_darts();
	m.edge_id(d, 0u);
	m.face_id(d, 0u);
	m.dart_level(d, m.current_level());

	// update max level if needed
	if(m.current_level() > m.maximum_level())
		m.maximum_level(m.current_level());
	return d;
}

}

#endif // MR_CMAP_INFOS_H
