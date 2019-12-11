#ifndef CGOGN_CORE_TYPES_CMAP_CPH_BASE_H
#define CGOGN_CORE_TYPES_CMAP_CPH_BASE_H

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/vector.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/cmap/cell.h>
#include <cgogn/core/functions/cmapbase_infos.h>

namespace cgogn{

template<typename MAP>
class CPH_Base{
public :
	using Self = CPH_Base;
	using AttributeContainer = AttributeContainerT<ChunkArray>;
	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	uint32 maximum_level_;

	/*!
	 * \brief Store the current number of darts per resolution level.
	 * This information is used to detect that a level of resolution
	 * becomes empty (contains no more dart) and that the maximum_level
	 * of the hierarchy should be decremented.
	 */
	std::vector<uint32> nb_darts_per_level_;
	std::shared_ptr<Attribute<uint32>> dart_level_;
	std::shared_ptr<MAP> base_map_;
public:

	inline CPH_Base() :
		maximum_level_(0u)
	{
		base_map_ = std::make_shared<MAP>();
		nb_darts_per_level_.reserve(32u);
		nb_darts_per_level_.push_back(0);
		dart_level_ = get_topology(*base_map_).template add_attribute<uint32>("dartLevel") ;
		if(dart_level_ == nullptr)
			dart_level_ = get_topology(*base_map_).template get_attribute<uint32>("dartLevel") ;
	}
	
	inline CPH_Base(std::shared_ptr<MAP>& m) :
		maximum_level_(0u),
		base_map_(m)
	{
		nb_darts_per_level_.reserve(32u);
		nb_darts_per_level_.push_back(0);
		dart_level_ = get_topology(*base_map_).template add_attribute<uint32>("dartLevel") ;
		if(dart_level_ == nullptr)
			dart_level_ = get_topology(*base_map_).template get_attribute<uint32>("dartLevel") ;
	}

	/***************************************************
	 *              LEVELS MANAGEMENT                  *
	 ***************************************************/

	inline uint32 maximum_level() const
	{
		return maximum_level_;
	}

	inline void maximum_level(uint32 l)
	{
		maximum_level_ = l;
	}

	inline uint32 dart_level(Dart d) const
	{
		return (*dart_level_)[d.index] ;
	}

	inline void dart_level(Dart d, uint32 l)
	{
		(*dart_level_)[d.index] = l ;
	}
	
	inline void change_dart_level(Dart d,uint32 l)
	{
		nb_darts_per_level_[dart_level(d)]--;
		nb_darts_per_level_[l]++;
		if(l>dart_level(d) && l>maximum_level()){
			maximum_level(l);
		}
		if(l<dart_level(d)){
			while(nb_darts_per_level_[maximum_level()] == 0u){
				--maximum_level_;
				nb_darts_per_level_.pop_back();
			}
		}
		(*dart_level_)[d.index] = l ;
	}
	
	inline const MAP& mesh() const {return *base_map_;}
	inline MAP& mesh(){return *base_map_;}
};

}


#endif // CGOGN_CORE_TYPES_CMAP_CPH_BASE_H
