#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/io/surface_import.h>

using namespace cgogn;

template <typename MESH>
void do_something(MESH& m)
{
	using Face = typename mesh_traits<MESH>::Face;

	Face f = add_face(m, 3);
	auto att = add_attribute<float64, Face>(m, "value");
	value<double>(m, att, f) = 2.0;
	foreach_cell(m, [&] (Face f)
	{
		std::cout << "face " << index_of(m, f) << " : " << value<float64>(m, att, f) << std::endl;
		return true;
	});
}



using Vec3 = Eigen::Vector3f;

int main()
{
	CMap2 map2;

	do_something(map2);

	std::cout << "nb darts: " << map2.nb_darts() << std::endl;
	std::cout << "vertex attributes:" << std::endl;
	for (auto ag : map2.attribute_containers_[CMap2::Vertex::ORBIT])
		std::cout << ag->name() << std::endl;

	io::import_OFF<Vec3>(map2, "/home/kraemer/Media/Data/surface/lowRes/cube_tri.off");

	std::cout << "nb darts: " << map2.nb_darts() << std::endl;
	std::cout << "vertex attributes:" << std::endl;
	for (auto ag : map2.attribute_containers_[CMap2::Vertex::ORBIT])
		std::cout << ag->name() << std::endl;

	auto position = get_attribute<Vec3, CMap2::Vertex>(map2, "position");

	foreach_cell(map2, [&] (CMap2::Vertex v) -> bool
	{
		const Vec3& vec = value<Vec3>(map2, position, v);
		std::cout << "vertex " << map2.embedding(v) << " : " << vec[0] << "," << vec[1] << "," << vec[2] << std::endl;
		return true;
	});

	remove_attribute<CMap2::Vertex>(map2, position);

	///////////////////////

	CMap1 map1;

	do_something(map1);

	std::cout << "nb darts: " << map1.nb_darts() << std::endl;
	for (auto ag : map1.attribute_containers_[CMap1::Face::ORBIT])
		std::cout << ag->name() << std::endl;

	auto att1 = get_attribute<float64, CMap2::Face>(map1, "value");

	foreach_cell(map1, [&] (CMap1::Face f) -> bool
	{
		std::cout << "face " << map1.embedding(f) << " : " << value<float64>(map1, att1, f) << std::endl;
		return true;
	});
}
