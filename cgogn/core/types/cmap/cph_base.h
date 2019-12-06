#ifndef CGOGN_CORE_TYPES_CMAP_CPH_BASE_H
#define CGOGN_CORE_TYPES_CMAP_CPH_BASE_H

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/vector.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/cmap/cell.h>

namespace cgogn{

class CPH_Base{
public :
	using Self = CPH_Base;
	using AttributeContainer = AttributeContainerT<ChunkArray>;
	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	uint32 current_level_;
	uint32 maximum_level_;

	/*!
	 * \brief Store the current number of darts per resolution level.
	 * This information is used to detect that a level of resolution
	 * becomes empty (contains no more dart) and that the maximum_level
	 * of the hierarchy should be decremented.
	 */
	std::vector<uint32> nb_darts_per_level_;
	std::shared_ptr<Attribute<uint32>> dart_level_;
public:

	inline CPH_Base(AttributeContainer& topology) :
		current_level_(0u),
		maximum_level_(0u)
	{
		nb_darts_per_level_.reserve(32u);
		nb_darts_per_level_.push_back(0);
		dart_level_ = topology.template add_attribute<uint32>("dartLevel") ;
		if(dart_level_ == nullptr)
			dart_level_ = topology.template get_attribute<uint32>("dartLevel") ;
	}

	/***************************************************
	 *              LEVELS MANAGEMENT                  *
	 ***************************************************/

	inline uint32 current_level() const
	{
		return current_level_;
	}

	inline void current_level(uint32 l)
	{
		current_level_ = l ;
	}

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

	inline void inc_current_level()
	{
		current_level_++;
	}

	inline void dec_current_level()
	{
		cgogn_message_assert(current_level_ > 0u, "dec_current_level : already at minimal resolution level");

		if (current_level_ == maximum_level_) {
			if (nb_darts_per_level_[current_level_] == 0u) {
				maximum_level_--;
				nb_darts_per_level_.pop_back();
			}
		}

		current_level_--;
	}

	inline void inc_nb_darts()
	{
		nb_darts_per_level_[current_level_]++;
	}
};

}


#endif // CGOGN_CORE_TYPES_CMAP_CPH_BASE_H
