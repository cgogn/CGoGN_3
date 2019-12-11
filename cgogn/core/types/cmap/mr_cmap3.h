#ifndef CGOGN_CORE_TYPES_CMAP_MR_CMAP3_H
#define CGOGN_CORE_TYPES_CMAP_MR_CMAP3_H

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/cmap/cph3.h>

namespace cgogn {

class MRCmap3{
public:
	using Self = MRCmap3;

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
	uint32 current_level_;
	
	static const bool is_mesh_view = true;
	
	std::shared_ptr<CPH3<CMap3>> base_cph_;

	MRCmap3():current_level_(0u){
		base_cph_ = std::make_shared<CPH3<CMap3>>();
	}
	MRCmap3(std::shared_ptr<CMap3> m):current_level_(0u){
		base_cph_ = std::make_shared<CPH3<CMap3>>(m);
	}
	MRCmap3(const std::shared_ptr<CPH3<CMap3>>& m):current_level_(0u){
		base_cph_ = m;
	}
	MRCmap3(const MRCmap3& mr2):MRCmap3(mr2.base_cph_){
		this->current_level(mr2.current_level());
	}
	MRCmap3(const MRCmap3& mr2,uint32 l):MRCmap3(mr2.base_cph_){
		this->current_level(l);
	}
	
	MRCmap3 get_new_level_view(uint32 l)const{
		MRCmap3 result(*this);
		result.current_level(l);
		return result;
	}
	
	inline uint32 current_level() const
	{
		return current_level_;
	}

	inline void current_level(uint32 l)
	{
		current_level_ = l ;
	}
	
	inline void inc_nb_darts()
	{
		base_cph_->nb_darts_per_level_[current_level_]++;
	}
	
	inline const CPH3<CMap3>& cph() const {return *base_cph_;}
	inline CPH3<CMap3>& cph(){return *base_cph_;}
	
	inline const CMap3& mesh() const {return cph().mesh();}
	inline CMap3& mesh(){return cph().mesh();}
};

}


#endif // CGOGN_CORE_TYPES_CMAP_MR_CMAP3_H
