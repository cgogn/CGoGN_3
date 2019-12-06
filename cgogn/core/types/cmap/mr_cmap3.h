#ifndef CGOGN_CORE_TYPES_CMAP_MR_CMAP3_H
#define CGOGN_CORE_TYPES_CMAP_MR_CMAP3_H

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/cmap/cph3.h>

namespace cgogn {

class MRCmap3 : public CPH3{
public:
	using Self = MRCmap3;
	using Inherit = CPH3;

	using Vertex = CMap3::Vertex;
	using Vertex2 = CMap3::Vertex2;
	using Edge = CMap3::Edge;
	using Edge2 = Cell<PHI2>;
	using Face = CMap3::Face;
	using Face2 = CMap3::Face2;
	using Volume = CMap3::Volume;
	using CC = CMap3::CC;

	using Cells = std::tuple<Vertex, Vertex2, Edge, Edge2, Face, Face2, Volume>;

	using AttributeContainer = AttributeContainerT<ChunkArray>;
	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	
	static const bool is_mesh_view = true;
	
	std::shared_ptr<CMap3> base_map_;

	MRCmap3():base_map_(),Inherit(base_map_->base_map_->topology_){}
	MRCmap3(std::shared_ptr<CMap3> m):base_map_(m),Inherit(m->base_map_->topology_){}

	MRCmap3(const MRCmap3& mr2):MRCmap3(mr2.base_map_){
		this->current_level(mr2.current_level());
	}
	MRCmap3(const MRCmap3& mr2,uint32 l):MRCmap3(mr2.base_map_){
		this->current_level(l);
	}
	
	inline const CMap3& mesh() const {return *base_map_;}
	inline CMap3& mesh(){return *base_map_;}
	
	MRCmap3 get_new_level_view(uint32 l)const{
		MRCmap3 result(*this);
		result.current_level(l);
		return result;
	}
};

}


#endif // CGOGN_CORE_TYPES_CMAP_MR_CMAP3_H
