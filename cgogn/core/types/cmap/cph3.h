#ifndef CPH3_H
#define CPH3_H

#include "cph2.h"

namespace cgogn
{

class CPH3 : public CPH2
{
public:

    using Self =  CPH3;
    using Inherit = CPH2;

    using AttributeContainer = AttributeContainerT<ChunkArray>;
    template <typename T>
    using Attribute = AttributeContainer::Attribute<T>;

protected:

    std::shared_ptr<Attribute<uint32>> face_id_;

public:
    CPH3(AttributeContainer& topology) : Inherit(topology)
    {
		face_id_ = topology.template add_attribute<uint32>("faceId");
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
};

} // namespace cgogn

#endif // CPH3_H
