/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#ifndef CGOGN_IO_SURFACE_OBJ_H_
#define CGOGN_IO_SURFACE_OBJ_H_

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/utils.h>

#include <fstream>
using namespace std::literals::string_literals;

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_OBJ(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	SurfaceImportData surface_data;

	std::ifstream fp(filename.c_str(), std::ios::in);

	std::string line, tag;
	line.reserve(512u);

	// read vertices position
	do
	{
		fp >> tag;
		if (tag == std::string("v"))
		{
			surface_data.nb_vertices_++;
			float64 x = read_double(fp, line);
			float64 y = read_double(fp, line);
			float64 z = read_double(fp, line);
			surface_data.vertex_position_.push_back({x, y, z});
		}
		getline_safe(fp, line); // flush line
	} while (!fp.eof());

	if (surface_data.nb_vertices_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no vertices." << std::endl;
		return false;
	}

	// rewind
	fp.clear();
	fp.seekg(0, std::ios::beg);

	// read faces (vertex indices)
	// uint32 nb_faces = 0;
	do
	{
		fp >> tag;
		getline_safe(fp, line);
		if (tag == std::string("f"))
		{
			surface_data.nb_faces_++;
			std::vector<uint32> indices;
			std::istringstream iss(line);
			std::string str;
			while (!iss.eof())
			{
				iss >> str;
				uint32 ind = 0;
				while ((ind < str.length()) && (str[ind] != '/'))
					ind++;
				if (ind > 0)
				{
					uint32 index;
					std::stringstream iss(str.substr(0, ind));
					iss >> index;
					indices.push_back(index - 1);
				}
			}
			surface_data.faces_nb_vertices_.push_back(indices.size());
			surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
													  indices.end());
		}
	} while (!fp.eof());

	if (surface_data.nb_faces_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no faces." << std::endl;
		return false;
	}

	import_surface_data(m, surface_data);

	return true;
}


template <typename MESH>
bool import_OBJ_tn(MESH& m_p, MESH& m_tc, MESH& m_n, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	Scoped_C_Locale loc;

	SurfaceImportData surface_data_p;
	SurfaceImportData2D surface_data_tc;
	SurfaceImportData surface_data_n;

	std::ifstream fp(filename.c_str(), std::ios::in);
	if (!fp.good())
	{
		std::cerr << "Error opening file " << filename << std::endl;
		return false;
	}
	fp.seekg(0, std::ios::end);
	uint64 sz = fp.tellg();
	fp.seekg(0, std::ios::beg);
	std::vector<char> buffer(sz + 1);
	fp.read(buffer.data(), sz);
	buffer[sz] = 0;
	std::string sbuffer(buffer.data());

	std::istringstream ss(sbuffer);

	std::string tag;
	std::string line;
	// read vertices position
	do
	{
		ss >> tag;
		if (tag == std::string("v"))
		{
			surface_data_p.nb_vertices_++;
			float64 x = read_double(ss, line);
			float64 y = read_double(ss, line);
			float64 z = read_double(ss, line);
			surface_data_p.vertex_position_.push_back({x, y, z});
		}

		if (tag == std::string("vt"))
		{
			surface_data_tc.nb_vertices_++;
			float64 x = read_double(ss, line);
			float64 y = read_double(ss, line);
			surface_data_tc.vertex_position_.push_back({x, y});
		}

		if (tag == std::string("vn"))
		{
			surface_data_n.nb_vertices_++;
			float64 x = read_double(ss, line);
			float64 y = read_double(ss, line);
			float64 z = read_double(ss, line);
			surface_data_n.vertex_position_.push_back({x, y, z});
		}
		getline_safe(ss, line); // flush line
	} while (!ss.eof());

	if (surface_data_p.nb_vertices_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no vertices." << std::endl;
		return false;
	}

	// rewind
	ss.clear();
	ss.seekg(0, std::ios::beg);

	// read faces (vertex indices)
	// uint32 nb_faces = 0;
	do
	{
		ss >> tag;
		getline_safe(ss, line);
		if (tag == "f"s)
		{
			std::vector<uint32> indices_p;
			std::vector<uint32> indices_tc;
			std::vector<uint32> indices_n;

			std::istringstream iss(line);
			while (!iss.eof())
			{
				std::string word_buf;
				iss >> word_buf;
				if (!word_buf.empty())
				{
					//					std::cout << "WORD: " << word_buf << std::endl;
					auto slash1 = word_buf.find('/');
					if (slash1 == std::string::npos)
					{
						uint32 index = std::atoi(word_buf.c_str());
						indices_p.push_back(index - 1);
					}
					else
					{
						auto slash2 = word_buf.find('/', slash1 + 1);
						if (slash2 == std::string::npos)
						{
							std::string str_ind = word_buf.substr(0, slash1);
							uint32 index = std::atoi(str_ind.c_str());
							indices_p.push_back(index - 1);
							str_ind = word_buf.substr(slash1 + 1, std::string::npos);
							index = std::atoi(str_ind.c_str());
							indices_tc.push_back(index - 1);
						}
						else
						{
							std::string str_ind = word_buf.substr(0, slash1);
							uint32 index = std::atoi(str_ind.c_str());
							indices_p.push_back(index - 1);
							if ((slash2 - slash1) > 1)
							{
								str_ind = word_buf.substr(slash1 + 1, slash2);
								index = std::atoi(str_ind.c_str());
								indices_tc.push_back(index - 1);
							}
							str_ind = word_buf.substr(slash2 + 1, std::string::npos);
							if (str_ind.size() > 0)
								index = std::atoi(str_ind.c_str());
							indices_n.push_back(index - 1);
						}
					}
				}
			}

			if (indices_p.size() > 0)
			{
				surface_data_p.nb_faces_++;
				surface_data_p.faces_nb_vertices_.push_back(indices_p.size());
				surface_data_p.faces_vertex_indices_.insert(surface_data_p.faces_vertex_indices_.end(),
															indices_p.begin(), indices_p.end());
			}


			if (indices_tc.size() > 0)
			{
				surface_data_tc.nb_faces_++;
				surface_data_tc.faces_nb_vertices_.push_back(indices_tc.size());
				surface_data_tc.faces_vertex_indices_.insert(surface_data_tc.faces_vertex_indices_.end(),
															 indices_tc.begin(), indices_tc.end());
			}

			if (indices_n.size() > 0)
			{
				surface_data_n.nb_faces_++;
				surface_data_n.faces_nb_vertices_.push_back(indices_n.size());
				surface_data_n.faces_vertex_indices_.insert(surface_data_n.faces_vertex_indices_.end(),
															indices_n.begin(), indices_n.end());
			}

		}
	} while (!ss.eof());

	if (surface_data_p.nb_faces_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no faces." << std::endl;
		return false;
	}

	import_surface_data(m_p, surface_data_p,true);
	import_surface_data(m_tc, surface_data_tc, false);
	import_surface_data(m_n, surface_data_n, false);
	bool mesh_tc_ok = (surface_data_tc.nb_faces_ > 0u);
	bool mesh_n_ok = (surface_data_n.nb_faces_ > 0u);
	
	if (mesh_tc_ok || mesh_n_ok)
	{
		uint32 nb_boundary_edges_tc = 0;
		uint32 nb_boundary_edges_n = 0;
		foreach_cell(m_p, [&](MESH::Edge e) -> bool {
			Dart d = e.dart_;
			Dart d2 = phi2(m_p,d);
			if (is_boundary(m_p, d2) || is_boundary(m_p, d))
				return true;
			if (d2 == d)
				return true;
			if (mesh_tc_ok)
			{
				if ((index_of(m_tc, MESH::Vertex(d2)) == index_of(m_tc, MESH::Vertex(phi1(m_tc, d)))) &&
					(index_of(m_tc, MESH::Vertex(d)) == index_of(m_tc, MESH::Vertex(phi1(m_tc, d2)))))
					phi2_sew(m_tc, d, d2);
				else
					nb_boundary_edges_tc++;
			}
			if (mesh_n_ok)
			{
				if ((index_of(m_n, MESH::Vertex(d2)) == index_of(m_n, MESH::Vertex(phi1(m_n, d)))) &&
					(index_of(m_n, MESH::Vertex(d)) == index_of(m_n, MESH::Vertex(phi1(m_n, d2)))))
					phi2_sew(m_n, d, d2);
				else
					nb_boundary_edges_n++;
			}
			return true;
		});
		uint32 nb_holes = close(m_tc);
		std::cout << nb_holes << " hole(s) have been closed in texture coord map" << std::endl;
		std::cout << nb_boundary_edges_tc << " edges in boundary of texture coord map" << std::endl;
		nb_holes = close(m_n);	
		std::cout << nb_holes << " hole(s) have been closed in normal map" << std::endl;
		std::cout << nb_boundary_edges_n << " edges in boundary of normal map" << std::endl;
	}


	return true;
}


template <typename MESH>
void export_OBJ(MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
				const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	// TODO
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_OBJ_H_
