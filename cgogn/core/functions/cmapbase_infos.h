#ifndef CGOGN_CORE_FUNCTIONS_CMAPBASE_INFOS_H
#define CGOGN_CORE_FUNCTIONS_CMAPBASE_INFOS_H

#include <cgogn/core/types/cmap/cmap_base.h>
#include <cgogn/core/types/cmap/mr_cmap3.h>

namespace cgogn {

//////////////
// CMapBase //
//////////////

inline bool is_boundary(const CMapBase& m, Dart d)
{
	return (*m.boundary_marker_)[d.index] != 0u;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH, typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline bool is_boundary(const MESH& m, Dart d)
{
	return is_boundary(m.mesh(),d);
}

//////////////
// CMapBase //
//////////////


inline void set_boundary(const CMapBase& m,Dart d, bool b)
{
	(*m.boundary_marker_)[d.index] = b ? 1u : 0u;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline void set_boundary(const MESH& m,Dart d, bool b)
{
	return set_boundary(m.mesh(),d,b);
}

//////////////
// CMapBase //
//////////////

inline Dart begin(const CMapBase& m){ return Dart(m.topology_.first_index()); }
inline Dart end(const CMapBase& m){ return Dart(m.topology_.last_index()); }
inline Dart next(const CMapBase& m,Dart d){ return Dart(m.topology_.next_index(d.index)); }

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline Dart begin(const MESH& m){ return begin(m.mesh()); }

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline Dart end(const MESH& m){ return end(m.mesh()); }

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline Dart next(const MESH& m,Dart d){ return next(m.mesh(),d); }

//////////////
// CMapBase //
//////////////

inline uint32 nb_darts(const CMapBase& m)
{
	return m.topology_.nb_elements();
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline uint32 nb_darts(const MESH& m)
{
	return nb_darts(m.mesh());
}

//////////////
// CMapBase //
//////////////

inline Dart add_dart(CMapBase& m)
{
	uint32 index = m.topology_.new_index();
	Dart d(index);
	for (auto rel : m.relations_)
		(*rel)[d.index] = d;
	for (auto emb : m.cells_indices_)
		if (emb)
			(*emb)[d.index] = INVALID_INDEX;
	return d;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline Dart add_dart(MESH& m)
{
	return add_dart(m.mesh());
}


}


#endif // CGOGN_CORE_FUNCTIONS_CMAPBASE_INFOS_H
