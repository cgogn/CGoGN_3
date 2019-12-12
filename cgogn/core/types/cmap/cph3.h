#ifndef CPH3_H
#define CPH3_H

#include "cph2.h"
#include <set>

namespace cgogn
{

template<typename MAP>
class CPH3 : public CPH2<MAP>
{
public:

    using Self =  CPH3;
    using Inherit = CPH2<MAP>;

    using AttributeContainer = AttributeContainerT<ChunkArray>;
    template <typename T>
    using Attribute = AttributeContainer::Attribute<T>;

protected:

    std::shared_ptr<Attribute<uint32>> face_id_;

public:
    CPH3() : Inherit()
    {
		face_id_ = get_topology(*this->base_map_).template add_attribute<uint32>("faceId");
    }
	
	inline CPH3(std::shared_ptr<MAP> m) : Inherit(m)
	{
		face_id_ = get_topology(*m).template add_attribute<uint32>("faceId") ;
		if(face_id_ == nullptr)
			face_id_ = get_topology(*m).template get_attribute<uint32>("faceId") ;
	}
    
    /***************************************************
     *             FACE ID MANAGEMENT                  *
     ***************************************************/
    
    inline uint32 face_id(Dart d) const
    {
		return (*face_id_)[d.index] ;
    }
    
    inline void face_id(Dart d, uint32 i)
    {
		(*face_id_)[d.index] = i ;
    }
	
	inline uint32 refinement_face_id(const std::vector<Dart>& cut_path) const
    {
		std::set<uint32> set_fid;
		for(Dart d:cut_path){
			set_fid.insert(face_id(d));
		}
		uint32 result = 0;
		while(set_fid.find(result) != set_fid.end()){
            ++result;
        }
		return result;
    }
};

} // namespace cgogn

#endif // CPH3_H
