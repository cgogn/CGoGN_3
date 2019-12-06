#ifndef CPH2_H
#define CPH2_H

#include "cph_base.h"

namespace cgogn 
{

class CPH2 : public CPH_Base
{
public:

    using Self = CPH2;
    using Inherit = CPH_Base;
    using AttributeContainer = AttributeContainerT<ChunkArray>;
    template <typename T>
    using Attribute = AttributeContainer::Attribute<T>;
	
protected:
    
    std::shared_ptr<Attribute<uint32>> edge_id_;
	
public:
	
    CPH2(AttributeContainer& topology) : Inherit(topology)
    {
		edge_id_ = topology.template add_attribute<uint32>("edgeId");
    }
    
    /***************************************************
    *             EDGE ID MANAGEMENT                  *
    ***************************************************/
    
    inline uint32 edge_id(Dart d) const
    {
		return (*edge_id_)[d.index] ;
    }
    
    inline void edge_id(Dart d, uint32 i)
    {
		(*edge_id_)[d.index] = i ;
    }
    
    inline uint32 refinement_edge_id(Dart d, Dart e) const
    {
		uint32 d_id = edge_id(d);
		uint32 e_id = edge_id(e);
		
		uint32 id = d_id + e_id;
		
		if (id == 0u)
			return 1u;
		else if (id == 1u)
			return 2u;
		else if (id == 2u)
		{
			if (d_id == e_id)
				return 0u;
			else
				return 1u;
		}
		// else if (id == 3)
		return 0u;
    }


};
}

#endif // CPH2_H
