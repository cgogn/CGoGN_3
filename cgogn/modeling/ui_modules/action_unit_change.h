/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_ACTION_UNIT_CHANGE_H_
#define CGOGN_ACTION_UNIT_CHANGE_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <thirdparty/rapidcsv/rapidcsv.h>

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#define PATH_OF_OPENFACE "../../../../OpenFace/build/bin/"
#define DEFAULT_PATH CGOGN_STR(CGOGN_DATA_PATH) "../../"

namespace fs = std::filesystem;

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

const Vec3 GREEN = Vec3(0, 128, 0);
const Vec3 BLUE = Vec3(0, 0, 255);

template <typename MESH>
class ActionUnitChange : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "ActionUnitChange can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

public:
	ActionUnitChange(const App& app)
		: ViewModule(app, "ActionUnitChange (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr), selected_vertex_position_(nullptr)
	{
	}
	~ActionUnitChange()
	{
	}

public:
	static bool ends_with(const std::string& str, const std::string& suffix)
	{
		return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}

	// call this function after you initialized the module
	void set_directory(std::string dirname)
	{
		directory_ = dirname;
	}

	// No signal system , call this function after you initialized the module
	void set_mesh(MESH& m, std::shared_ptr<Attribute<Vec3>> vertex_position)
	{
		selected_mesh_ = &m;
		selected_vertex_position_ = vertex_position;
	}

	// No signal system , call this function after you initialized the module
	void set_view(View& v)
	{
		selected_view_ = &v;
	}

	// Put inside a vector all the files with an extension ext
	void set_all_paths(std::string root, std::string ext, std::vector<std::string>& paths)
	{
		for (auto& p : fs::recursive_directory_iterator(root))
		{
			if (p.path().extension() == ext)
			{
				paths.push_back(p.path().string());
				std::cout << p.path().string() << std::endl;
			}
		}
		std::sort(paths.begin(), paths.end());
	}

	// Create a new attribute or get an attribute and fill it with data from another attribute
	void set_attribute(MESH& m, Attribute<Vec3>* to_set, std::string attribute_name, float weight)
	{
		std::shared_ptr<Attribute<Vec3>> attribute_to_change =
			cgogn::get_or_add_attribute<Vec3, Vertex>(m, attribute_name.c_str());
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, attribute_to_change, v) = value<Vec3>(m, to_set, v) * weight;
			return true;
		});
	}

	void change_to_selected_au(MESH& m, Attribute<Vec3>* au_position)
	{
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");
		Attribute<Vec3>* vertex_pos_value = vertex_position.get();
		Vec3 tmp = Vec3(0,0,0);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, vertex_pos_value, v) = value<Vec3>(m, au_position, v);
			tmp = value<Vec3>(m, au_position, v);
			return true;
		});
		mesh_provider_->emit_attribute_changed(m, vertex_pos_value);
	}

	// Function for taking a screenshot of the face to send to OpenFace for analysis
	void take_screenshot(int num, std::string dir_of_name)
	{
		std::ostringstream name;
		name << DEFAULT_PATH << "CGoGN_3/build/stage/bin/";
		name << "Screenshot_";
		for (int j = 0; j < 4 - std::to_string(num).size(); j++)
		{
			name << "0";
		}
		name << num;
		name << ".jpg";

		std::ostringstream dirname;
		dirname << DEFAULT_PATH << "OpenFace/samples/" << dir_of_name << "/";

		selected_view_->save_screenshot_name(name.str());
		fs::path sourceFile = name.str().c_str();
		fs::path targetParent = dirname.str().c_str();
		if (!fs::is_directory(targetParent) || !fs::exists(targetParent))
			fs::create_directory(targetParent);

		fs::copy(sourceFile, targetParent, fs::copy_options::overwrite_existing);
	}

	// This function creates all the differents AUs that have been found with set_all_paths and create for each of them
	// an attribute DO NOT USE LOAD_SURFACE_FROM_FILE since it creates a new mesh and causes problems with the signal
	// system
	void setup_mesh_attributes()
	{
		for (auto path : path_aus_)
		{
			std::cout << path.substr(path.size() - 8, path.size() - (path.size() - 8) - 4) << std::endl;
			std::ifstream fp(path.c_str(), std::ios::in);
			if (!fp.good())
			{
				std::cerr << "Error opening file " << path.c_str() << std::endl;
				return;
			}
			std::shared_ptr<Attribute<Vec3>> au_pos = cgogn::add_attribute<Vec3, Vertex>(
				*selected_mesh_, path.substr(path.size() - 8, path.size() - (path.size() - 8) - 4));
			pos_aus_.push_back(au_pos);
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
			std::vector<Vec3> vec_pos;
			// read vertices position
			do
			{
				ss >> tag;
				if (tag == std::string("v"))
				{
					float64 x = cgogn::io::read_double(ss, line);
					float64 y = cgogn::io::read_double(ss, line);
					float64 z = cgogn::io::read_double(ss, line);
					Vec3 temp = Vec3(x, y, z);
					vec_pos.push_back(temp);
				}
			} while (!ss.eof());

			int incr = 0;
			Vec3 point_norm;
			cgogn::foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
				point_norm = vec_pos[index_of(*selected_mesh_, v)];
				value<Vec3>(*selected_mesh_, au_pos, v) = point_norm;
				incr++;
				return true;
			});

			geometry::rescale(*au_pos, 1);
			mesh_provider_->emit_attribute_changed(*selected_mesh_, au_pos.get());
		}
	}

	//
	void setup_csv_matrix(int incr, std::shared_ptr<Attribute<Vec3>> au, float weight)
	{
		blending(*selected_mesh_, {au}, {weight});
		if (incr == pos_aus_.size())
		{
			std::ostringstream command;
			command << PATH_OF_OPENFACE << "csv_script_matrix.sh" << " " << DEFAULT_PATH << "OpenFace/samples/Matrix/"
					<< " " << directory_ << "CSV/" << " " << DEFAULT_PATH << "OpenFace/build/bin/" << " "
					<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/rewrite_csv.py";
			if (system(command.str().c_str()) == 0)
			{
				if (system("rm -rf *.jpg") == 0)
					std::cout << "Erasing screenshots" << std::endl;
				else
					std::cout << "Command invalid" << std::endl;

				std::ostringstream matrix_csv_path;
				matrix_csv_path << directory_ << "CSV/Matrix.csv";
				std::string path = matrix_csv_path.str();
				csv_parser(path, ',', csv_weights_detected_, csv_weights_confirm_, vector_OF_rest_csv_);
				set_matrix_jacob();
			}
			else
				std::cout << "Command invalid" << std::endl;
		}
	}

	// Define AU00 as vector at rest
	void setup_vector_at_rest()
	{
		std::ostringstream command;
		command << PATH_OF_OPENFACE << "csv_script_matrix.sh" << " " << DEFAULT_PATH << "OpenFace/samples/Matrix/"
				<< " " << directory_ << "CSV/" << " " << DEFAULT_PATH << "OpenFace/build/bin/";

		take_screenshot(0, "Matrix");
		std::cout << command.str().c_str() << std::endl;
		if (system(command.str().c_str()) == 0)
		{
			std::ostringstream matrix_csv_path;
			matrix_csv_path << directory_ << "CSV/Matrix.csv";
			std::string path = matrix_csv_path.str();
			csv_parser(path, ',', csv_weights_detected_, csv_weights_confirm_, vector_OF_rest_csv_);
			vector_OF_rest_cgogn_ = vector_OF_rest_csv_;
			std::cout << "Vector of rest cgogn " << std::endl << vector_OF_rest_cgogn_ << std::endl;
			matrix_jacob.resize(csv_weights_confirm_.cols(), csv_weights_confirm_.cols());
		}
		else
			std::cout << "Command invalid" << std::endl;
	}

	// Calculate the jacobian matrix for changing space between OpenFace space and CGoGN space
	void set_matrix_jacob()
	{
		matrix_jacob.resize(csv_weights_detected_.rows(), csv_weights_detected_.cols());

		for (int i = 0; i < csv_weights_detected_.rows(); i++)
		{
			for (int j = 0; j < csv_weights_detected_.cols(); j++)
			{
				float element_jacob =
					(csv_weights_detected_(i, j) - vector_OF_rest_cgogn_(j)) / weight_for_jacob_matrix;
				matrix_jacob(i, j) = ((abs(element_jacob) > epsilon) ? element_jacob : 0);
			}
		}

		std::cout << "Jacobian Matrix : " << std::endl << matrix_jacob.format(OctaveFmt) << std::endl;

		std::cout << "Vector from face at rest : " << std::endl << vector_OF_rest_cgogn_ << std::endl;
	}

	// Apply the jacobian Matrix to selected CSV
	// Can use the confirmation of weights
	void apply_matrix_csv(bool confirm)
	{
		Eigen::VectorXd tmp;
		tmp.resize(csv_weights_detected_.cols());
		for (int i = 0; i < csv_weights_detected_.rows(); i++)
		{
			for (int j = 0; j < csv_weights_detected_.cols(); j++)
			{
				if (csv_weights_confirm_(i, j) == 0 && confirm)
					csv_weights_detected_(i, j) = 0.;
				else
				{
					csv_weights_detected_(i, j) = csv_weights_detected_(i, j) + vector_OF_rest_cgogn_(j);
				}
			}
			tmp = matrix_jacob * csv_weights_detected_.row(i).transpose();
			csv_weights_detected_.row(i) = tmp.transpose();
		}

		std::cout << "Matrix of detection : " << std::endl << csv_weights_detected_.format(OctaveFmt) << std::endl;
	}

	// Function to compare original csv and csv passed inside the OpenFace
	// Do not use since there's more frame than normal in the new csv
	void compare_matrices(Eigen::MatrixXd& csv_weights_detected, Eigen::MatrixXd& compare_weights_detected)
	{
		for (int i = 0; i < compare_weights_detected.rows(); i++)
		{
			for (int j = 0; j < compare_weights_detected.cols(); j++)
			{
				compare_weights_detected(i, j) = (compare_weights_detected(i, j) - vector_OF_rest_csv_(j));
				compare_weights_detected(i, j) = compare_weights_detected(i, j) - csv_weights_detected(i, j);
			}
		}
		std::cout << "Matrice de la comparaison entre le csv et le nouveau csv : " << std::endl
				  << compare_weights_detected.format(OctaveFmt) << std::endl;
	}

	// Creates either the jacobian matrix or a test matrix to verify that the jacobian isn't false
	void create_matrix_test_and_jacob(Eigen::MatrixXd& matrix, int& incr, bool resize, bool jacob)
	{
		take_screenshot(0, "test");
		std::ostringstream command;
		command << PATH_OF_OPENFACE << "csv_script_matrix.sh" << " " << DEFAULT_PATH << "OpenFace/samples/test/" << " "
				<< directory_ << "CSV/" << " " << DEFAULT_PATH << "OpenFace/build/bin/" << " "
				<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/rewrite_csv.py";
		if (system(command.str().c_str()) == 0)
		{
			std::ostringstream path;
			path << directory_ << "CSV/test.csv";
			std::string string_path = path.str();
			Eigen::MatrixXd weights_detected_;
			Eigen::MatrixXd weights_confirm_;
			Eigen::VectorXd vec;
			csv_parser(string_path, ',', weights_detected_, weights_confirm_, vec);
			if (resize)
				matrix.resize(weights_detected_.cols(), weights_detected_.cols());

			if (jacob)
			{
				for (int j = 0; j < weights_detected_.cols(); j++)
				{
					float element_jacob =
						(weights_detected_(0, j) - vector_OF_rest_cgogn_(j)) / weight_for_jacob_matrix;
					std::cout << element_jacob << std::endl;
					matrix_jacob(incr - 2, j) = ((abs(element_jacob) > epsilon) ? element_jacob : 0);
				}
			}
			else
			{
				for (int j = 0; j < weights_detected_.cols(); j++)
				{
					float element_jacob = (weights_detected_(0, j) - vector_OF_rest_cgogn_(j));
					matrix(incr - 2, j) = element_jacob;
				}
			}
		}
	}

	// Function that will take the different screenshots and do the different blendings to calculate test matrices
	// Will write inside a named file the results of tests
	void loop_openface(Eigen::MatrixXd& matrix, int& incr, bool& resize, float weight, bool& start, bool confidence,
					   std::string name)
	{
		if ((incr == pos_aus_.size()))
		{
			create_matrix_test_and_jacob(matrix, incr, resize, false);
			test_matrix(matrix_jacob, matrix, weight, name, vector_confidence_lower_bound,
						vector_confidence_upper_bound, confidence);
			incr = 1;
			start = false;
		}
		if ((incr < pos_aus_.size()))
		{
			if (incr == 1)
			{
				setup_csv_matrix(incr, pos_aus_[incr], weight);
				incr++;
				resize = true;
			}
			else
			{
				setup_csv_matrix(incr, pos_aus_[incr], weight);
				create_matrix_test_and_jacob(matrix, incr, resize, false);
				if (resize)
					resize = false;
				incr++;
			}
		}
	}

	// Function to calculate the lower bound and upper bound of the jacobian matrix since OpenFace is a neural network
	// it's necessary to do to validate the results
	void calculate_confidence_interval(Eigen::MatrixXd& matrix1, Eigen::MatrixXd& matrix2, Eigen::MatrixXd& jacob,
									   float weight_used1, float weight_used2,
									   Eigen::VectorXd& confidence_lower_bound_vector,
									   Eigen::VectorXd& confidence_upper_bound_vector, std::string name)
	{
		confidence_lower_bound_vector.resize(jacob.cols());
		confidence_upper_bound_vector.resize(jacob.cols());

		Eigen::MatrixXd confidence_lower_bound_jacob;
		Eigen::MatrixXd confidence_upper_bound_jacob;

		Eigen::MatrixXd tmp = jacob;

		confidence_lower_bound_jacob.resize(jacob.rows(), jacob.cols());
		confidence_upper_bound_jacob.resize(jacob.rows(), jacob.cols());

		confidence_lower_bound_jacob.setZero(jacob.rows(), jacob.rows());
		confidence_upper_bound_jacob.setZero(jacob.rows(), jacob.rows());

		std::ofstream outputFile(name);

		if (outputFile.is_open())
		{
			for (int i = 0; i < matrix1.rows(); i++)
			{
				float element_jacob = matrix1(i, i) / weight_used1;
				confidence_lower_bound_jacob(i, i) = ((abs(element_jacob) > epsilon) ? element_jacob : 0);
				element_jacob = matrix2(i, i) / weight_used2;
				confidence_upper_bound_jacob(i, i) = ((abs(element_jacob) > epsilon) ? element_jacob : 0);

				if ((tmp(i, i) < confidence_lower_bound_jacob(i, i)) &&
					(tmp(i, i) < confidence_upper_bound_jacob(i, i)))
				{
					if (confidence_lower_bound_jacob(i, i) < confidence_upper_bound_jacob(i, i))
					{
						float new_value_for_slop = confidence_lower_bound_jacob(i, i);
						confidence_lower_bound_jacob(i, i) = tmp(i, i);
						tmp(i, i) = new_value_for_slop;
					}
					else
					{
						float new_value_for_slop = confidence_upper_bound_jacob(i, i);
						confidence_upper_bound_jacob(i, i) = tmp(i, i);
						tmp(i, i) = new_value_for_slop;
					}
				}

				if ((tmp(i, i) > confidence_lower_bound_jacob(i, i)) &&
					(tmp(i, i) > confidence_upper_bound_jacob(i, i)))
				{
					if (confidence_lower_bound_jacob(i, i) > confidence_upper_bound_jacob(i, i))
					{
						float new_value_for_slop = confidence_lower_bound_jacob(i, i);
						confidence_lower_bound_jacob(i, i) = tmp(i, i);
						tmp(i, i) = new_value_for_slop;
					}
					else
					{
						float new_value_for_slop = confidence_upper_bound_jacob(i, i);
						confidence_upper_bound_jacob(i, i) = tmp(i, i);
						tmp(i, i) = new_value_for_slop;
					}
				}

				if (confidence_lower_bound_jacob(i, i) > confidence_upper_bound_jacob(i, i))
				{
					float inversion = confidence_upper_bound_jacob(i, i);
					confidence_upper_bound_jacob(i, i) = confidence_lower_bound_jacob(i, i);
					confidence_lower_bound_jacob(i, i) = inversion;
				}
			}

			jacob = tmp;

			for (int i = 0; i < matrix1.rows(); i++)
			{
				if (jacob(i, i) != 0)
				{
					confidence_lower_bound_vector(i) = confidence_lower_bound_jacob(i, i) - jacob(i, i);
					confidence_upper_bound_vector(i) = confidence_upper_bound_jacob(i, i) - jacob(i, i);
				}
				else
				{
					confidence_lower_bound_vector(i) = confidence_lower_bound_jacob(i, i);
					confidence_upper_bound_vector(i) = confidence_upper_bound_jacob(i, i);
				}
			}
			outputFile << std::endl;
			outputFile << "tmp matrix : " << std::endl << tmp.format(OctaveFmt) << std::endl;
			outputFile << "matrix jacobian : " << std::endl << jacob.format(OctaveFmt) << std::endl;

			outputFile << "Lower Bound confidence : " << std::endl;
			outputFile << "Confidence matrix for weight " << weight_used1 << " : " << std::endl
					   << matrix1.format(OctaveFmt) << std::endl;
			outputFile << "Confidence matrix jacobian for weight " << weight_used1 << " : " << std::endl
					   << confidence_lower_bound_jacob.format(OctaveFmt) << std::endl;
			outputFile << "Confidence Vector for weight " << weight_used1 << " : " << std::endl
					   << confidence_lower_bound_vector.format(VectorFmt) << std::endl;

			outputFile << "Upper Bound confidence : " << std::endl;
			outputFile << "Confidence matrix for weight " << weight_used2 << " : " << std::endl
					   << matrix2.format(OctaveFmt) << std::endl;
			outputFile << "Confidence matrix jacobian for weight " << weight_used2 << " : " << std::endl
					   << confidence_upper_bound_jacob.format(OctaveFmt) << std::endl;
			outputFile << "Confidence Vector for weight " << weight_used2 << " : " << std::endl
					   << confidence_upper_bound_vector.format(VectorFmt) << std::endl;
		}
		outputFile.close();
	}

	// This function will write in the file name all the unit tests for a matrice while comparing the results with the
	// jacobian if needed it can not use the confidence interval Will write for the diagonal if the AU passed the test
	void test_matrix(Eigen::MatrixXd& jacobian, Eigen::MatrixXd& matrix_test, float weight_matrix, std::string name,
					 Eigen::VectorXd confidence_lower_bound, Eigen::VectorXd confidence_upper_bound,
					 bool use_confidence)
	{
		Eigen::MatrixXd res;
		int nb_pass = 0;
		int nb_failed = 0;
		int nb_total_pass = 0;
		int nb_total_failed = 0;
		int nb_test_pass_specific_AU = 0;
		int nb_test_fail_specific_AU = 0;
		res.resize(jacobian.rows(), jacobian.cols());
		float value_with_lower_bound = 0.;
		float value_with_upper_bound = 0.;

		std::ofstream outputFile(name);
		if (outputFile.is_open())
		{
			if (use_confidence)
			{
				outputFile << "Using confidence interval : " << std::endl;
				outputFile << confidence_lower_bound.format(OctaveFmt) << std::endl;
				outputFile << confidence_upper_bound.format(OctaveFmt) << std::endl;
			}
			else
				outputFile << "Without confidence interval : " << std::endl;

			for (int i = 0; i < res.rows(); i++)
			{
				nb_pass = 0;
				nb_failed = 0;

				for (int j = 0; j < res.cols(); j++)
				{
					if (use_confidence)
					{
						res(i, j) = (jacobian(i, j) * weight_matrix) - matrix_test(i, j);
						value_with_lower_bound =
							((jacobian(i, j) + confidence_lower_bound(j)) * weight_matrix) - matrix_test(i, j);
						value_with_upper_bound =
							((jacobian(i, j) + confidence_upper_bound(j)) * weight_matrix) - matrix_test(i, j);

						if (((value_with_lower_bound <= res(i, j) && value_with_upper_bound >= res(i, j))))
						{
							std::cout << "Where : " << i << " " << j << " : value : " << res(i, j)
									  << " lower bound : " << confidence_lower_bound(i, j) << " "
									  << " upper_bound : " << confidence_upper_bound(i, j) << std::endl;
						}

						if (abs(res(i, j)) > 0.1 &&
							(value_with_lower_bound > res(i, j) || value_with_upper_bound < res(i, j)))
							nb_failed++;
						else
							nb_pass++;

						if ((i == j) && (abs(res(i, j)) > 0.1) &&
							(value_with_lower_bound > res(i, j) || value_with_upper_bound < res(i, j)))
						{
							outputFile << "Test failed for : " << pos_aus_[i + 1]->name().c_str() << std::endl;
							nb_test_fail_specific_AU++;
						}
						else if ((i == j) && ((abs(res(i, j)) < 0.1) || (value_with_lower_bound <= res(i, j) &&
																		 value_with_upper_bound >= res(i, j))))
						{
							outputFile << "Test passed for : " << pos_aus_[i + 1]->name().c_str() << std::endl;
							nb_test_pass_specific_AU++;
						}
					}
					else
					{
						res(i, j) = (jacobian(i, j) * weight_matrix) - matrix_test(i, j);

						if (abs(res(i, j)) > 0.1 &&
							(value_with_lower_bound > res(i, j) || value_with_upper_bound < res(i, j)))
							nb_failed++;
						else
							nb_pass++;

						if (i == j && abs(res(i, j)) > 0.1)
						{
							outputFile << "Test failed for : " << pos_aus_[i + 1]->name().c_str() << std::endl;
							nb_test_fail_specific_AU++;
						}
						else if (i == j && abs(res(i, j)) < 0.1)
						{
							outputFile << "Test passed for : " << pos_aus_[i + 1]->name().c_str() << std::endl;
							nb_test_pass_specific_AU++;
						}
					}
				}

				nb_total_failed += nb_failed;
				nb_total_pass += nb_pass;
				outputFile << "Number of tests passed  : " << nb_pass << std::endl;
				outputFile << "Number of tests failed  : " << nb_failed << std::endl << std::endl;
				outputFile << "Ratio of tests passed  : " << (float(nb_pass) / float((nb_failed + nb_pass))) * 100
						   << "%" << std::endl
						   << std::endl;
			}
			outputFile << std::endl;
			outputFile << "Res Matrix  : " << std::endl << res.format(OctaveFmt) << std::endl;
			outputFile << "Jacobian Matrix  : " << std::endl << jacobian.format(OctaveFmt) << std::endl;
			outputFile << "Matrix testing  : " << std::endl << matrix_test.format(OctaveFmt) << std::endl;

			outputFile << std::endl;
			outputFile << "Total Number of tests : " << nb_total_failed + nb_total_pass << std::endl;
			outputFile << "Total Number of tests passed : " << nb_total_pass << std::endl;
			outputFile << "Total Number of tests failed : " << nb_total_failed << std::endl;
			outputFile << "Ratio of total tests passed : "
					   << (float(nb_total_pass) / float((nb_total_failed + nb_total_pass))) * 100 << "%" << std::endl;

			outputFile << std::endl;
			outputFile << "Total Number of tests for the diagonal of the matrix : "
					   << nb_test_fail_specific_AU + nb_test_pass_specific_AU << std::endl;
			outputFile << "Total Number of tests passed for the diagonal of the matrix : " << nb_test_pass_specific_AU
					   << std::endl;
			outputFile << "Total Number of tests failed for the diagonal of the matrix : " << nb_test_fail_specific_AU
					   << std::endl;
			outputFile << "Ratio of tests passed for the diagonal of the matrix : "
					   << (float(nb_test_pass_specific_AU) /
						   float((nb_test_fail_specific_AU + nb_test_pass_specific_AU))) *
							  100
					   << "%" << std::endl;

			outputFile << std::endl;
			outputFile << "Determinant of jacobian matrix : " << jacobian.determinant() << std::endl;

			outputFile << "Inverse" << std::endl << matrix_jacob.inverse().format(OctaveFmt) << std::endl;
		}
		std::cout << std::endl;
		std::cout << "Text has been written to the file." << std::endl;
		outputFile.close();
	}

	// This function will write in the file name all the unit tests for a line while comparing the results with the
	// jacobian TO DO , change to a vector of weights for combinations of AUs with differents weights Will write for
	// each AU in the face if the AU passed the test
	void test_line(Eigen::MatrixXd& jacobian, Eigen::VectorXd& line_test, Eigen::VectorXd confidence_lower_bound,
				   Eigen::VectorXd confidence_upper_bound, std::vector<int> au_in, float weight_au, std::string name)
	{
		Eigen::VectorXd res;
		Eigen::VectorXd res_lower_value;
		Eigen::VectorXd res_upper_value;
		int nb_pass = 0;
		int nb_failed = 0;
		int nb_test_pass_specific_AU = 0;
		int nb_test_fail_specific_AU = au_in.size();
		res.resize(jacobian.cols());
		res_upper_value.resize(jacobian.cols());
		res_lower_value.resize(jacobian.cols());
		line_test = line_test - vector_OF_rest_cgogn_;

		std::ofstream outputFile(name);
		if (outputFile.is_open())
		{
			for (int j = 0; j < res.rows(); j++)
			{
				float tmp = 0.;
				float value_with_lower_bound = 0.;
				float value_with_upper_bound = 0.;
				for (int i = 0; i < au_in.size(); i++)
				{
					tmp += jacobian(au_in[i], j) * weight_au;
					value_with_lower_bound += (jacobian(au_in[i], j) + confidence_lower_bound(j)) * weight_au;
					value_with_upper_bound += (jacobian(au_in[i], j) + confidence_upper_bound(j)) * weight_au;
				}

				res(j) = tmp - line_test(j);
				res_lower_value(j) = value_with_lower_bound - line_test(j);
				res_upper_value(j) = value_with_upper_bound - line_test(j);
				if (abs(res(j)) > 0.1 && (res_lower_value(j) > res(j) || res_upper_value(j) < res(j)))
					nb_failed++;
				else
					nb_pass++;

				if ((std::find(std::begin(au_in), std::end(au_in), j)) != std::end(au_in) &&
					((abs(res(j)) < 0.1) || (res_lower_value(j) <= res(j) && res_upper_value(j) >= res(j))))
				{
					outputFile << "Test passed for : " << pos_aus_[j]->name().c_str() << std::endl;
					nb_test_pass_specific_AU++;
					nb_test_fail_specific_AU--;
				}
			}
			outputFile << name << std::endl;
			outputFile << "Number of tests passed  : " << nb_pass << std::endl;
			outputFile << "Number of tests failed  : " << nb_failed << std::endl << std::endl;
			outputFile << "Ratio of tests passed  : " << (float(nb_pass) / float((nb_failed + nb_pass))) * 100 << "%"
					   << std::endl
					   << std::endl;

			outputFile << std::endl;
			outputFile << "Res Line  : " << std::endl << res.transpose().format(OctaveFmt) << std::endl;
			outputFile << "Res Upper Value Line  : " << std::endl
					   << res_upper_value.transpose().format(OctaveFmt) << std::endl;
			outputFile << "Res Lower Value Line  : " << std::endl
					   << res_lower_value.transpose().format(OctaveFmt) << std::endl;
			outputFile << "Line test : " << std::endl << line_test.transpose().format(OctaveFmt) << std::endl;
			outputFile << std::endl;

			std::cout << "Res Upper Value Line  : " << std::endl
					  << res_upper_value.transpose().format(OctaveFmt) << std::endl;
			std::cout << "Res Lower Value Line  : " << std::endl
					  << res_lower_value.transpose().format(OctaveFmt) << std::endl;

			outputFile << "Confidence vector upper value : " << std::endl
					   << confidence_upper_bound.transpose().format(OctaveFmt) << std::endl;
			outputFile << "Confidence vector lower value  : " << std::endl
					   << confidence_lower_bound.transpose().format(OctaveFmt) << std::endl;
			outputFile << "Jacobian Matrix  : " << std::endl << jacobian.format(OctaveFmt) << std::endl;
			outputFile << std::endl;

			outputFile << std::endl;
			outputFile << "Total Number of tests for the AUs used in the line : "
					   << nb_test_fail_specific_AU + nb_test_pass_specific_AU << std::endl;
			outputFile << "Total Number of tests passed for the AUs used in the line : " << nb_test_pass_specific_AU
					   << std::endl;
			outputFile << "Total Number of tests failed for the AUs used in the line : " << nb_test_fail_specific_AU
					   << std::endl;
			outputFile << "Ratio of tests passed for the AUs used in the line : "
					   << (float(nb_test_pass_specific_AU) /
						   float((nb_test_fail_specific_AU + nb_test_pass_specific_AU))) *
							  100
					   << "%" << std::endl;
		}

		std::ostringstream command;
		command << "rm -rf " << directory_ << "CSV/" << name << "*";

		if (system(command.str().c_str()) == 0)
		{
			std::cout << "Removing CSV file" << std::endl;
		}

		std::cout << std::endl;
		std::cout << "Text has been written to the file " << name << std::endl;
		outputFile.close();
	}

	// Parse a csv using a filename and the separator of the csv
	void csv_parser(std::string& filename, char separator, Eigen::MatrixXd& csv_weights_detected,
					Eigen::MatrixXd& csv_weights_confirm, Eigen::VectorXd& vector_OF_rest_csv)
	{
		rapidcsv::Document doc(filename, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams(separator, true));
		std::ofstream outputFile("test.txt"); // Open/create a file named "test.txt" for writing
		std::vector<std::string> csv_columns_name = doc.GetColumnNames();
		std::vector<float> tempData;
		int nb_elem = 0;
		csv_.clear();
		for (int i = 0; i < csv_columns_name.size(); i++)
		{
			if (strcmp(csv_columns_name[i].c_str(), "AU28_c") != 0)
			{
				std::cout << csv_columns_name[i].c_str() << std::endl;
				tempData = doc.GetColumn<float>(csv_columns_name[i]);
				for (int j = 0; j < tempData.size(); j++)
				{
					if (outputFile.is_open())
					{
						outputFile << tempData[j];
						outputFile << ";";
					}
					else
					{
						std::cout << "Failed to create the file." << std::endl;
					}
				}
				outputFile << tempData.size();
				outputFile << "\n";
				csv_.emplace(csv_columns_name[i], tempData);
				tempData.clear();
			}
		}
		std::cout << "Text has been written to the file." << std::endl;
		outputFile.close();

		int i = 0;
		int nb_columns = 0;
		std::map<std::string, std::vector<float>>::iterator iter = csv_.begin();
		nb_elem = iter->second.size();
		for (auto& it : csv_)
		{
			if (ends_with(it.first, "_r") || ends_with(it.first, "_c"))
			{
				nb_columns++;
			}
		}

		int incr = 0;
		int incr2 = 0;

		if (nb_elem > 1)
		{
			nb_elem--;
		}

		std::cout << "NB_columns : " << nb_columns / 2 << std::endl;
		std::cout << "NB_elems : " << nb_elem << std::endl;
		csv_weights_detected.resize(nb_elem, nb_columns / 2);
		csv_weights_confirm.resize(nb_elem, nb_columns / 2);
		vector_OF_rest_csv.resize(nb_columns / 2);
		vector_OF_rest_cgogn_.resize(nb_columns / 2);

		for (auto& it : csv_)
		{
			if (ends_with(it.first, "_r"))
			{
				vector_OF_rest_csv(incr) = it.second[0];
				if (it.second.size() <= 1)
					csv_weights_detected(0, incr) = it.second[0];

				for (int j = 1; j < it.second.size(); j++)
					csv_weights_detected(j - 1, incr) = it.second[j];
				incr++;
			}
			if (ends_with(it.first, "_c"))
			{
				if (it.second.size() <= 1)
					csv_weights_confirm(0, incr) = it.second[0];
				for (int j = 1; j < it.second.size(); j++)
					csv_weights_confirm(j - 1, incr2) = it.second[j];
				incr2++;
			}
		}
		std::cout << "Vector of rest csv : " << std::endl << vector_OF_rest_csv_ << std::endl;
		std::cout << "Matrix of detection : " << std::endl << csv_weights_detected.format(OctaveFmt) << std::endl;
		std::cout << "Matrix of confirmation : " << std::endl << csv_weights_confirm.format(OctaveFmt) << std::endl;
	}

	// write jacobian matrix inside a file to not have to redo all the loop to calculate it
	void write_jacob_to_file(std::string name, Eigen::MatrixXd& jacobian, Eigen::VectorXd confidence_lower_bound,
							 Eigen::VectorXd confidence_upper_bound)
	{
		std::ofstream outputFile(name);
		if (outputFile.is_open())
		{
			outputFile << std::setprecision(2);
			outputFile << jacobian.rows() << " " << jacobian.cols() << std::endl;
			outputFile << jacobian << std::endl;
			outputFile << confidence_lower_bound.rows() << std::endl;
			outputFile << confidence_lower_bound << std::endl;
			outputFile << confidence_upper_bound << std::endl;
			outputFile << vector_OF_rest_cgogn_ << std::endl;
		}
		outputFile.close();
	}

	// read the file with the jacobian matrix in it
	// if it works , will disable the loop to create the jacobian matrix
	bool read_jacob_from_file(std::string name)
	{
		bool jacob_read = false;
		std::string line;
		std::ifstream inputFile;
		inputFile.open(name);
		if (inputFile.is_open())
		{
			int jacob_nb_rows = 0;
			int jacob_nb_cols = 0;
			inputFile >> jacob_nb_rows;
			inputFile >> jacob_nb_cols;
			matrix_jacob.resize(jacob_nb_rows, jacob_nb_cols);

			for (int i = 0; i < jacob_nb_rows; i++)
			{
				for (int j = 0; j < jacob_nb_cols; j++)
				{
					inputFile >> matrix_jacob(i, j);
				}
			}

			int vector_size = 0;
			inputFile >> vector_size;
			vector_confidence_lower_bound.resize(vector_size);
			vector_confidence_upper_bound.resize(vector_size);
			vector_OF_rest_cgogn_.resize(vector_size);

			for (int i = 0; i < vector_size; i++)
			{
				inputFile >> vector_confidence_lower_bound(i);
			}

			for (int i = 0; i < vector_size; i++)
			{
				inputFile >> vector_confidence_upper_bound(i);
			}

			for (int i = 0; i < vector_size; i++)
			{
				inputFile >> vector_OF_rest_cgogn_(i);
			}
			jacob_read = true;
		}
		inputFile.close();

		std::cout << matrix_jacob.rows() << " " << matrix_jacob.cols() << std::endl;
		std::cout << matrix_jacob.format(OctaveFmt) << std::endl;
		std::cout << vector_confidence_lower_bound.rows() << std::endl;
		std::cout << vector_confidence_lower_bound.format(OctaveFmt) << std::endl;
		std::cout << vector_confidence_upper_bound.format(OctaveFmt) << std::endl;
		std::cout << vector_OF_rest_cgogn_.format(OctaveFmt) << std::endl;

		return jacob_read;
	}

	// Set the color of the points to blue
	// Use this function before changing
	void set_to_blue(MESH& m)
	{
		std::shared_ptr<Attribute<Vec3>> color_change = cgogn::get_attribute<Vec3, Vertex>(m, "color");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, color_change, v) = BLUE;
			return true;
		});
		mesh_provider_->emit_attribute_changed(m, color_change.get());
	}

	// Change the color of the point to green
	// Used for optimisation (WIP)
	void highlight_difference(MESH& m, std::shared_ptr<Attribute<Vec3>> au_position)
	{
		std::shared_ptr<Attribute<Vec3>> color_change = cgogn::get_attribute<Vec3, Vertex>(m, "color");
		std::shared_ptr<Attribute<Vec3>> pos_au_repos = cgogn::get_attribute<Vec3, Vertex>(m, "AU00");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (value<Vec3>(m, au_position, v) != value<Vec3>(m, pos_au_repos, v))
				value<Vec3>(m, color_change, v) = GREEN;
			else
				value<Vec3>(m, color_change, v) = BLUE;
			return true;
		});
		mesh_provider_->emit_attribute_changed(m, color_change.get());
	}

	// Blending function
	// Currently it's a sum of vectors of each different attributes that will be blent
	void blending(MESH& m, std::vector<std::shared_ptr<Attribute<Vec3>>> attributes_to_blend, std::vector<float> weight_list)
	{
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");
		std::shared_ptr<Attribute<Vec3>> color = cgogn::get_attribute<Vec3, Vertex>(m, "color");
		std::shared_ptr<Attribute<Vec3>> repos_position = cgogn::get_attribute<Vec3, Vertex>(m, "AU00");
		Attribute<Vec3>* new_vertex_pos_value = vertex_position.get();
		Vec3 diff_distance_repos = Vec3(0, 0, 0);
		float epsilon = 0.0001;

		// need to check if parallel_foreach_cell messes with the calculations 
		foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, vertex_position, v) = value<Vec3>(m, repos_position, v);
			Vec3 result = Vec3(0, 0, 0);
			float nb_au_influence = 0.;

			for (int i = 0; i < attributes_to_blend.size(); i++)
			{	
				diff_distance_repos = value<Vec3>(m, attributes_to_blend[i], v) - value<Vec3>(m, repos_position, v);
				
				if(abs(diff_distance_repos[0]) > 0. || abs(diff_distance_repos[1]) > 0. || abs(diff_distance_repos[2]) > 0.){
					nb_au_influence++;
				}
				result += diff_distance_repos * weight_list[i];
			}

			if (nb_au_influence != 0.)
				result = result / nb_au_influence;

			result[0] = (abs(result[0]) > epsilon) ? result[0] : 0. ; 
			result[1] = (abs(result[1]) > epsilon) ? result[1] : 0. ;
			result[2] = (abs(result[2]) > epsilon) ? result[2] : 0. ;
			
			value<Vec3>(m, vertex_position, v) += result ;
			return true;
		});

		mesh_provider_->emit_attribute_changed(m, new_vertex_pos_value);
	}

	// Function called each frame after clicking on Apply CSV
	// Setup the differents AUs and weights to calculate the frames
	void blending_csv(std::vector<float> weights, int incr, float poids_frame)
	{
		std::vector<std::shared_ptr<Attribute<Vec3>>> attributes_csv;
		std::vector<float> weight_list;

		for (auto& it : csv_)
		{
			if (ends_with(it.first, "_r"))
			{
				attributes_csv.push_back(
					cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, it.first.substr(0, it.first.size() - 2)));
			}
		}
		for (int i = 1; i < pos_aus_.size(); i++)
		{
			if (incr + 1 < csv_weights_detected_.rows())
			{
				weight_list.push_back((csv_weights_detected_(incr, i - 1) * (1 - poids_frame)) +
						 (csv_weights_detected_(incr + 1, i - 1) * poids_frame));
			}
			else
				weight_list.push_back((csv_weights_detected_(incr - 1, i - 1) * (1 - poids_frame)) +
						 (csv_weights_detected_(incr, i - 1) * poids_frame));
		}
		blending(*selected_mesh_, attributes_csv, weight_list);
	}

	// Compute the distance between points in the starting configuration and the end configuration for the interpolation
	// algorithm
	void set_distance(MESH& m, Attribute<Vec3>* blendshape_start, Attribute<Vec3>* blendshape_target,
					  float weight_start, float weight_target)
	{
		std::shared_ptr<Attribute<Vec3>> distance = cgogn::get_or_add_attribute<Vec3, Vertex>(m, "distance");
		std::shared_ptr<Attribute<Vec3>> repos_position = cgogn::get_attribute<Vec3, Vertex>(m, "AU00");
		Attribute<Vec3>* distance_value = distance.get();
		Vec3 diff_distance_repos = Vec3(0, 0, 0);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			diff_distance_repos =
				((value<Vec3>(m, blendshape_target, v) - value<Vec3>(m, repos_position, v)) * weight_target) -
				((value<Vec3>(m, blendshape_start, v) - value<Vec3>(m, repos_position, v)) * weight_start);
			value<Vec3>(m, distance_value, v) = diff_distance_repos;
			return true;
		});
	}

	// Interpolation function with a step
	// Call set_distance before this function
	// Switch attribute to position_interpolation to watch the interpolation
	void interpolation(MESH& m, float pas)
	{
		std::shared_ptr<Attribute<Vec3>> distance = cgogn::get_or_add_attribute<Vec3, Vertex>(m, "distance");
		std::shared_ptr<Attribute<Vec3>> position_interpolation =
			cgogn::get_or_add_attribute<Vec3, Vertex>(m, "position_interpolation");
		Attribute<Vec3>* interpolation_value = position_interpolation.get();

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, position_interpolation, v) = value<Vec3>(m, distance, v) * pas;
			return true;
		});
		mesh_provider_->emit_attribute_changed(m, interpolation_value);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		set_all_paths(directory_, ".obj", path_aus_);
		set_all_paths(directory_, ".csv", path_csv_);
		jacob_read = read_jacob_from_file("jacob.txt");
	}

	void left_panel() override
	{
		ImGui::InputText("FPS", std::to_string(ui::App::fps()).data(), std::to_string(ui::App::fps()).size());

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			selected_vertex_position_.reset();
			// selected_vertex_normal_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			static float weight = 0.;
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			// imgui_combo_attribute<Vertex, Vec3>(
			// 	*selected_mesh_, selected_vertex_normal_, "Normal",
			// 	[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_normal_ = attribute; });

			ImGui::Separator();

			if (selected_vertex_position_)
			{
				change_to_selected_au(*selected_mesh_, selected_vertex_position_.get());

				static const char* current_item_mesh = NULL;
				if (ImGui::BeginCombo(
						"Meshes to blend",
						current_item_mesh))
				{

					for (int n = 0; n < pos_aus_.size(); n++)
					{
						bool is_selected = (current_item_mesh ==
											pos_aus_[n]->name().c_str()); 
						if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected))
						{
							current_item_mesh = pos_aus_[n]->name().c_str();
							if (!(std::find(std::begin(attribute_to_blend_), std::end(attribute_to_blend_),
											pos_aus_[n]) != std::end(attribute_to_blend_)))
							{
								attribute_to_blend_.push_back(pos_aus_[n]);
								weights.push_back(1.);
							}
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					}
					ImGui::EndCombo();
				}

				static bool start = false;
				static int i = 0;
				static int nb_au = 1;
				if (attribute_to_blend_.size() != 0)
				{
					for (int i = 0; i < attribute_to_blend_.size(); i++)
					{
						ImGui::SliderFloat(attribute_to_blend_[i]->name().c_str(), &weights[i], -1.0, 5.0);
					}

					if (ImGui::Button("Blend"))
					{
						std::ostringstream new_attribute_name;
		
						for (int i = 0; i < attribute_to_blend_.size(); i++)
						{
							new_attribute_name << attribute_to_blend_[i]->name().c_str() << 'w' << weights[i] << '+';
						}
						blending(*selected_mesh_, attribute_to_blend_, weights);

						std::shared_ptr<Attribute<Vec3>> new_attribute = add_attribute<Vec3, Vertex>(
							*selected_mesh_,
							new_attribute_name.str().substr(0, new_attribute_name.str().size() - 1).c_str());

						// start = true;
						attribute_to_blend_.clear();
						weights.clear();
					}
				}

				if (ImGui::Button("Blend progressif"))
				{
					highlight_difference(*selected_mesh_, pos_aus_[1]);
					start = true;
				}

				if (ImGui::Button("Clear"))
				{
					std::shared_ptr<Attribute<Vec3>> repos_position =
						cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, "AU00");
					blending(*selected_mesh_, {repos_position}, {weight});
					start = false;
					weight = -1.;
					attribute_to_blend_.clear();
				}

				if (ImGui::Button("Stop"))
				{
					start = false;
				}

				if (start)
				{
					if (weight < 5.)
					{
						blending(*selected_mesh_, {pos_aus_[nb_au]}, {weight});
						weight += 0.003;
						i++;
					}
					else
					{
						nb_au++;
						if (nb_au >= pos_aus_.size())
						{
							start = false;
							std::ostringstream command;
							command << PATH_OF_OPENFACE << "script.sh" << " " << DEFAULT_PATH << "OpenFace/samples/"
									<< " " << directory_ << "CSV_VIDEO/" << " " << DEFAULT_PATH << "OpenFace/build/bin/"
									<< " "
									<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/parser.py";
							if (system(command.str().c_str()) == 0)
								std::cout << "Command succesfully executed" << std::endl;
							else
								std::cout << "Command impossible to execute" << std::endl;
						}
						else
						{
							weight = 0.;
							blending(*selected_mesh_, {pos_aus_[nb_au]}, {weight});
							i = 0;
						}
					}
				}

				ImGui::Separator();

				static const char* current_item_csv = NULL;
				if (ImGui::BeginCombo("Load CSV", current_item_csv))
				{
					for (int n = 0; n < path_csv_.size(); n++)
					{
						bool is_selected = (current_item_csv == path_csv_[n].c_str());
						if (ImGui::Selectable(path_csv_[n].substr(directory_.size(), path_csv_[n].size()).c_str(),
											  is_selected))
						{
							current_item_csv = path_csv_[n].c_str();
							csv_.clear();
							timestamp_csv_.clear();
							csv_parser(path_csv_[n], ',', csv_weights_detected_, csv_weights_confirm_,
									   vector_OF_rest_csv_);
							apply_matrix_csv(false);
							i = 0;
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					}
					ImGui::EndCombo();
				}

				static int incr = 0;
				static int nb_screen = 0;
				static float poids_frame = 1.;
				if (current_item_csv != NULL)
				{
					if (ImGui::Button("Apply CSV"))
					{
						std::map<std::string, std::vector<float>>::iterator iter = csv_.begin();
						count_timer_csv = iter->second.size();
						for (auto& it : csv_)
						{
							if (it.first == "timestamp")
								timestamp_csv_ = it.second;
						}
						for (int i = 0; i < timestamp_csv_.size(); i++)
						{
							std::cout << timestamp_csv_[i] << std::endl;
						}

						time_start = ui::App::frame_time_;
						incr = 0.;
						poids_frame = 1.;
					}

					if (ImGui::Button("Stop CSV"))
					{
						incr = count_timer_csv + 1;
					}
				}

				if (incr < count_timer_csv)
				{
					blending_csv(weights, incr, poids_frame);
					timer = ui::App::frame_time_ - time_start;

					while ((timer > timestamp_csv_[incr]) && (incr < count_timer_csv))
					{
						incr++;
					}
					if (incr - 1 > 0)
					{
						poids_frame =
							(timer - timestamp_csv_[incr - 1]) / (timestamp_csv_[incr] - timestamp_csv_[incr - 1]);
					}
					// if (nb_screen > 0)
					// {
					// 	take_screenshot(nb_screen - 1, "CSV");
					// }
					std::cout << "timer : " << timer << std::endl;
					std::cout << "poids_frame : " << poids_frame << std::endl;
					std::cout << "timestamp_csv : " << timestamp_csv_[incr] << std::endl;
					nb_screen++;
				}

				ImGui::Separator();

				static const char* current_item_start = NULL;
				if (ImGui::BeginCombo(
						"Interpolation shape start",
						current_item_start)) // The second parameter is the label previewed before opening the combo.
				{

					for (int n = 0; n < pos_aus_.size(); n++)
					{
						bool is_selected = (current_item_start ==
											pos_aus_[n]->name().c_str()); // You can store your selection however you
																		  // want, outside or inside your objects
						if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected))
						{
							current_item_start = pos_aus_[n]->name().c_str();
							au_start = pos_aus_[n].get();
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus(); // You may set the initial focus when opening the combo
														  // (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}

				if (current_item_start != NULL)
				{
					ImGui::SliderFloat(au_start->name().c_str(), &weight_start, 0.0, 5.0);
				}

				static const char* current_item_target = NULL;
				if (ImGui::BeginCombo(
						"Interpolation shape target",
						current_item_target)) // The second parameter is the label previewed before opening the combo.
				{

					for (int n = 0; n < pos_aus_.size(); n++)
					{
						bool is_selected = (current_item_target ==
											pos_aus_[n]->name().c_str()); // You can store your selection however you
																		  // want, outside or inside your objects
						if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected))
						{
							current_item_target = pos_aus_[n]->name().c_str();
							au_target = pos_aus_[n].get();
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus(); // You may set the initial focus when opening the combo
														  // (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}

				if (current_item_target != nullptr)
				{
					ImGui::SliderFloat(au_target->name().c_str(), &weight_target, 0.0, 5.0);
				}

				ImGui::SliderInt("Nb_frames", &nb_frames, 120, 600);

				if (ImGui::Button("Start Interpolation"))
				{
					if ((current_item_start != NULL) && (current_item_target != NULL))
					{
						set_attribute(*selected_mesh_, au_start, "position_interpolation", weight_start);
						set_distance(*selected_mesh_, au_start, au_target, weight_start, weight_target);
						count_interpolation = nb_frames;
					}
				}

				if (count_interpolation != 0)
				{
					interpolation(*selected_mesh_, 1. / nb_frames);
					count_interpolation--;
				}

				if (start)
				{
					take_screenshot(i, pos_aus_[nb_au]->name());
				}

				static int matrix_incr = 0;
				static int matrix_incr1 = 1;
				static int matrix_incr2 = 1;
				static bool resize = false;
				static bool test = false;
				static bool confidence_1 = false;
				static bool confidence_2 = false;
				static bool multi_au_test = false;
				const float weight_confidence1 = weight_for_jacob_matrix - 0.25;
				const float weight_confidence2 = weight_for_jacob_matrix + 0.5;
				static Eigen::MatrixXd matrix_confidence_1;
				static Eigen::MatrixXd matrix_confidence_2;

				if (!jacob_read)
				{
					if (matrix_incr == pos_aus_.size())
					{
		
						create_matrix_test_and_jacob(matrix_jacob, matrix_incr, false, true);
						matrix_incr++;
						confidence_1 = true;
						std::cout << "Jacobian Matrix : " << std::endl << matrix_jacob.format(OctaveFmt) << std::endl;
						std::cout << "Vector from face at rest : " << std::endl << vector_OF_rest_cgogn_ << std::endl;
					}
					if (matrix_incr < pos_aus_.size())
					{
						if (matrix_incr == 0)
						{
			
							setup_vector_at_rest();
							matrix_incr++;
			
							setup_csv_matrix(matrix_incr, pos_aus_[matrix_incr], weight_for_jacob_matrix);
						}
						else
						{
			
							setup_csv_matrix(matrix_incr, pos_aus_[matrix_incr], weight_for_jacob_matrix);
							if (matrix_incr >= 2)
								create_matrix_test_and_jacob(matrix_jacob, matrix_incr, false, true);
							matrix_incr++;
						}
					}

					if (confidence_1)
					{
						loop_openface(matrix_confidence_1, matrix_incr1, resize, weight_confidence1, confidence_1,
									  false, "./test_first_confidence_matrix.txt");
						if (!confidence_1)
						{
							confidence_2 = true;
						}
					}

					if (confidence_2)
					{
						loop_openface(matrix_confidence_2, matrix_incr2, resize, weight_confidence2, confidence_2,
									  false, "./test_second_confidence_matrix.txt");
						if (!confidence_2)
						{
							calculate_confidence_interval(matrix_confidence_1, matrix_confidence_2, matrix_jacob,
														  weight_confidence1, weight_confidence2,
														  vector_confidence_lower_bound, vector_confidence_upper_bound,
														  "./results_confidence.txt");
			
							std::cout << "Jacobian Matrix " << std::endl << matrix_jacob.format(OctaveFmt) << std::endl;
							write_jacob_to_file("jacob.txt", matrix_jacob, vector_confidence_lower_bound,
												vector_confidence_upper_bound);
						}
					}
				}

				static int nb_test = 5;
				static int incr_nb_test = 0.;
				static int nb_au_to_blend = 2;
				static std::vector<int> au_used;
				static std::string name_mix_au;

				ImGui::Separator();

				ImGui::SliderInt("Nb random AUs to blend", &nb_au_to_blend, 1, 10);
				ImGui::SliderInt("Number of tests", &nb_test, 1, 10);

				if (ImGui::Button("Multiples AUs test"))
				{
					multi_au_test = true;
					incr_nb_test = nb_test;
				}

				if (multi_au_test && incr_nb_test < nb_test)
				{
					take_screenshot(0, name_mix_au);
					std::ostringstream command;
					command << PATH_OF_OPENFACE << "csv_script_matrix.sh" << " " << DEFAULT_PATH << "OpenFace/samples/"
							<< name_mix_au << "/" << " " << directory_ << "CSV/" << " " << DEFAULT_PATH
							<< "OpenFace/build/bin/" << " "
							<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/rewrite_csv.py";
					if (system(command.str().c_str()) == 0)
					{
						std::ostringstream path;
						path << directory_ << "CSV/" << name_mix_au << ".csv";
						std::string string_path = path.str();
						Eigen::MatrixXd compare_weights_detected_;
						Eigen::MatrixXd compare_weights_confirm_;
						Eigen::VectorXd vec_compare_rest;
						csv_parser(string_path, ',', compare_weights_detected_, compare_weights_confirm_,
								   vec_compare_rest);
						test_line(matrix_jacob, vec_compare_rest, vector_confidence_lower_bound,
								  vector_confidence_upper_bound, au_used, 1., name_mix_au);
						au_used.clear();
						name_mix_au.clear();
					}

					if (incr_nb_test == 0)
						multi_au_test = false;
				}
				if (multi_au_test && incr_nb_test <= nb_test)
				{
					std::random_device rd;
					std::mt19937 gen(rd());
					std::uniform_int_distribution<> distr(1, pos_aus_.size() - 1);
					std::ostringstream name;
					std::vector<std::shared_ptr<Attribute<Vec3>>> attributes_au_used;
					std::vector<float> weights_used;
	
					multi_au_test = true;
					int nb_rand = 0;

					for (int i = 0; i < nb_au_to_blend; i++)
					{
						nb_rand = distr(gen);
						while (std::find(std::begin(au_used), std::end(au_used), nb_rand) != std::end(au_used))
						{
							nb_rand = distr(gen);
						}
						au_used.push_back(nb_rand);
						attributes_au_used.push_back(pos_aus_[nb_rand]);
						weights_used.push_back(1.);
						name << pos_aus_[nb_rand]->name().c_str() << "+";
						//highlight_difference(*selected_mesh_, pos_aus_[nb_rand]);
					}
					blending(*selected_mesh_, attributes_au_used, weights_used);
					name_mix_au = name.str().substr(0, name.str().size() - 1);
					incr_nb_test--;
				}

				ImGui::Separator();

				static int test_incr = 1;
				static float weight_for_test = 1.;

				ImGui::SliderFloat("Weight for test", &weight_for_test, 0.0, 5.0);

				if (ImGui::Button("Test jacobian"))
				{
					test = true;
				}

				if (test)
					loop_openface(matrix_test_unique_AU, test_incr, resize, weight_for_test, test, true,
								  "./results_test_jacobian.txt");

				static bool test_jacob = false;
				if (ImGui::Button("Test jacobian matrix"))
				{
					test_jacob = true;
				}

				if (test_jacob)
				{
					static int increment_matrix = 1;
					static Eigen::MatrixXd matrix_res_for_jacobian;

					if (increment_matrix == pos_aus_.size())
					{
		
						create_matrix_test_and_jacob(matrix_res_for_jacobian, increment_matrix, false, false);
						test_matrix(matrix_jacob, matrix_res_for_jacobian, 1, "./results_test_jacobian.txt",
									vector_confidence_lower_bound, vector_confidence_upper_bound, true);
						test_jacob = false;
						increment_matrix = 0;
						std::cout << "Matrix with each line of jacobian used : " << std::endl
								  << matrix_res_for_jacobian.format(OctaveFmt) << std::endl;
						std::cout << "Jacobian: " << std::endl << matrix_jacob.format(OctaveFmt) << std::endl;

						matrix_jacob = matrix_jacob.inverse();
					}
					if (increment_matrix < pos_aus_.size())
					{
						if (increment_matrix == 1)
						{
							matrix_jacob = matrix_jacob.inverse();
			
							// for (int i = 1; i < pos_aus_.size(); i++)
							// {
							//
							// 	blending(*selected_mesh_, pos_aus_[i], matrix_jacob(increment_matrix - 1, i - 1));
							// }
							//blending(*selected_mesh_, pos_aus_, matrix_jacob.row(increment_matrix - 1));
							increment_matrix++;
						}
						else
						{
			
							// for (int i = 1; i < pos_aus_.size(); i++)
							// {
							//
							// 	blending(*selected_mesh_, pos_aus_[i], matrix_jacob(increment_matrix - 1, i - 1));
							// }
							//blending(*selected_mesh_, pos_aus_, matrix_jacob.row(increment_matrix - 1));
							if (increment_matrix >= 2)
							{
								create_matrix_test_and_jacob(matrix_res_for_jacobian, increment_matrix, resize, false);
							}
							if (resize)
							{
								resize = false;
							}
							increment_matrix++;
						}
					}
				}
			}
		}
	}

private:
	MESH* selected_mesh_;
	View* selected_view_;
	MeshProvider<MESH>* mesh_provider_;

	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> pos_aus_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> attribute_to_blend_;

	Attribute<Vec3>* au_target;
	Attribute<Vec3>* au_start;

	std::vector<std::string> path_aus_;
	std::vector<std::string> path_csv_;
	std::string directory_;
	std::map<std::string, std::vector<float>> csv_;
	Eigen::MatrixXd csv_weights_detected_;
	Eigen::MatrixXd csv_weights_confirm_;

	Eigen::MatrixXd matrix_jacob;
	Eigen::MatrixXd matrix_test_unique_AU;

	Eigen::IOFormat OctaveFmt = Eigen::IOFormat(2, 0, ", ", ";\n", "", "", "[", "]");
	Eigen::IOFormat VectorFmt = Eigen::IOFormat(4, 0, ", ", ";\n", "", "", "[", "]");

	std::vector<float> weights;
	std::vector<float> timestamp_csv_;

	Eigen::VectorXd vector_OF_rest_cgogn_;
	Eigen::VectorXd vector_OF_rest_csv_;

	Eigen::VectorXd vector_confidence_lower_bound;
	Eigen::VectorXd vector_confidence_upper_bound;

	int nb_frames = 1;
	float64 timer = 0.;
	float64 time_start = 0.;
	int count_interpolation = 0;
	int count_timer_csv = 0;
	float weight_start = 1.;
	float weight_target = 1.;
	float weight_for_jacob_matrix = 1.5;
	float epsilon = 0.01;
	bool jacob_read = false;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_ACTION_UNIT_CHANGE_H_
