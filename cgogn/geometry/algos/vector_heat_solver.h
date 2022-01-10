#ifndef CGOGN_MODULE_VECTOR_HEAT_SOLVER_H_
#define CGOGN_MODULE_VECTOR_HEAT_SOLVER_H_

#include <Eigen/Sparse>
#include <Eigen/Dense>


#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/core/functions/attributes.h>

namespace cgogn
{

namespace geometry
{


template <typename MESH>
class VectorHeatSolver{

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using Edge = typename mesh_traits<MESH>::Edge;

private:
    Eigen::SparseMatrix<Scalar> massMatrix;
    //Eigen::SparseMatrix Lc;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>> scalarHeatSolver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<Scalar>>> vectorHeatSolver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>> poissonSolver;
    double shortTime;
    uint32 nb_vertices;

    std::shared_ptr<Attribute<uint32>> vertex_index;
    MESH* _mesh;

public:

VectorHeatSolver(MESH& mesh, Attribute<Vec3>* position, double time_multiplier = 1.0){

        vertex_index = get_or_add_attribute<uint32, Vertex>(mesh, "__vector_heat_index");
        this->_mesh = &mesh;

        nb_vertices = 0;
        foreach_cell(mesh, [&](Vertex v) -> bool {
            value<uint32>(mesh, vertex_index, v) = nb_vertices++;
            return true;
        });

		auto area = get_or_add_attribute<Scalar, Vertex>(mesh, "__area");
		geometry::compute_area<Vertex>(mesh, position, area.get());


		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc = geometry::cotan_operator_matrix(mesh, vertex_index.get(), position);


		Eigen::VectorXd vec_areas = Eigen::VectorXd(nb_vertices);

		foreach_cell(mesh, [&](Vertex v) -> bool {
			vec_areas(value<uint32>(mesh, vertex_index.get(), v)) =
				value<Scalar>(mesh, area, v);
			return true;
		});

		massMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor>(vec_areas.asDiagonal());

		// compute short time
		double t = 0;
		int nb_edges = 0;
		foreach_cell(mesh, [&](Edge e) -> bool {
			t += geometry::length(mesh, e, position);
			nb_edges++;
			return true;
		});

		t /= nb_edges;
		t *= t;
		t *= time_multiplier;
        shortTime = t;

		Eigen::SparseMatrix<Scalar> heatOp = massMatrix + t*Lc;
		scalarHeatSolver.compute(heatOp);

		// Vertex Connection Laplacian
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lconn = topo_laplacian_matrix(mesh, vertex_index.get());

		Eigen::SparseMatrix<std::complex<Scalar>> vectorOp = massMatrix.cast<std::complex<Scalar>>() + t*Lconn;


		//vectorHeatSolver.compute(vectorOp); TODO

		poissonSolver.compute(Lc);
}

~VectorHeatSolver(){
    remove_attribute<Vertex>(*_mesh, vertex_index);
}

};

}

}

#endif
