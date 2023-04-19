
#include <cgogn/io/surface/off.h>

#include <cgogn/geometry/types/vector_traits.h>

using namespace cgogn;

template <typename MESH>
void do_something(MESH& m)
{
	using Face = typename mesh_traits<MESH>::Face;

	Face f = add_face(m, 3);
	auto att = add_attribute<float64, Face>(m, "value");
	value<double>(m, att, f) = 2.0;
	foreach_cell(m, [&](Face f) {
		std::cout << "face " << index_of(m, f) << " : " << value<float64>(m, att, f) << std::endl;
		return true;
	});
}

using Vec3 = geometry::Vec3;

int main()
{
	CMap2 map2;

	do_something(map2);

	std::cout << "nb darts: " << nb_darts(map2) << std::endl;
	std::cout << "vertex attributes:" << std::endl;
	for (auto ag : map2.attribute_containers_[CMap2::Vertex::ORBIT])
		std::cout << ag->name() << std::endl;
	std::cout << "face attributes:" << std::endl;
	for (auto ag : map2.attribute_containers_[CMap2::Face::ORBIT])
		std::cout << ag->name() << std::endl;

	io::import_OFF(map2, "/home/kraemer/Media/Data/surface/lowRes/cube_tri.off");

	std::cout << "nb darts: " << nb_darts(map2) << std::endl;
	std::cout << "vertex attributes:" << std::endl;
	for (auto ag : map2.attribute_containers_[CMap2::Vertex::ORBIT])
		std::cout << ag->name() << std::endl;

	auto bla = add_attribute<uint32, CMap2::Vertex>(map2, "bla");
	auto bli = get_attribute<uint32, CMap2::Vertex>(map2, "bla");

	foreach_cell(map2, [&](CMap2::Vertex v) -> bool {
		uint32 i = value<uint32>(map2, bli, v);
		uint32 j = value<uint32>(map2, bla, v);
		std::cout << "vertex " << index_of(map2, v) << " : " << i << "," << j << std::endl;
		return true;
	});

	auto position = get_attribute<Vec3, CMap2::Vertex>(map2, "position");
	if (!position)
		std::cout << "position not valid" << std::endl;

	foreach_cell(map2, [&](CMap2::Vertex v) -> bool {
		const Vec3& vec = value<Vec3>(map2, position, v);
		std::cout << "vertex " << index_of(map2, v) << " : " << vec[0] << "," << vec[1] << "," << vec[2] << std::endl;
		return true;
	});

	remove_attribute<CMap2::Vertex>(map2, position);

	std::cout << "vertex attributes:" << std::endl;
	for (auto ag : map2.attribute_containers_[CMap2::Vertex::ORBIT])
		std::cout << ag->name() << std::endl;

	///////////////////////

	CMap1 map1;

	do_something(map1);

	std::cout << "nb darts: " << nb_darts(map1) << std::endl;
	std::cout << "face attributes:" << std::endl;
	for (auto ag : map1.attribute_containers_[CMap1::Face::ORBIT])
		std::cout << ag->name() << std::endl;

	auto att1 = get_attribute<float64, CMap2::Face>(map1, "value");

	foreach_cell(map1, [&](CMap1::Face f) -> bool {
		std::cout << "face " << index_of(map1, f) << " : " << value<float64>(map1, att1, f) << std::endl;
		return true;
	});

	//	remove_attribute<CMap1::Face>(map1, att1);
}
