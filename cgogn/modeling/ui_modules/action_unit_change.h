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
#include <cgogn/geometry/algos/picking.h>
#include <thirdparty/rapidcsv/rapidcsv.h>

#include <cgogn/modeling/algos/blending.h>
#include <cgogn/rendering/shape_drawer.h>

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#define DEFAULT_PATH CGOGN_STR(CGOGN_DATA_PATH) "../../"

namespace fs = std::filesystem;

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec3f;
using geometry::Vec4;
using GLMat4 = Eigen::Matrix4f;

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

	std::vector<std::shared_ptr<Attribute<Vec3>>> pos_aus_;

public:
	static bool ends_with(const std::string& str, const std::string& suffix)
	{
		return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}

	// 
	void exec_mode(bool use_exec_mode , std::string name_csv_to_use, std::string path_video)
	{
		use_exec_mode_ = use_exec_mode;
		exec_mode_ = use_exec_mode;
		exec_csv_name = name_csv_to_use;
		path_video_ = path_video;
	}

	// call this function after you initialized the module
	void set_directory(std::string dirname)
	{
		directory_ = dirname;
	}

	// call this function after you initialized the module
	void set_pathOpenface(std::string path_openface)
	{
		path_openface_ = path_openface;
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

	// No signal system , call this function after you initialized the module
	void set_position_attr_name(std::string name)
	{
		pos_attr_name = name;
	}

	// Put inside a vector all the files with an extension ext
	void set_all_paths(std::string root, std::string ext, std::vector<std::string>& paths)
	{
		for (auto& p : fs::recursive_directory_iterator(root))
		{
			if (p.path().extension() == ext)
				paths.push_back(p.path().string());
		}
		std::sort(paths.begin(), paths.end());
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
		dirname << path_openface_ << "samples/" << dir_of_name << "/";

		selected_view_->save_screenshot_name(name.str());
		fs::path sourceFile = name.str().c_str();
		fs::path targetParent = dirname.str().c_str();
		if (!fs::is_directory(targetParent) || !fs::exists(targetParent))
			fs::create_directory(targetParent);

		fs::copy(sourceFile, targetParent, fs::copy_options::overwrite_existing);
	}

	// Generate scripts needed to create the jacobian matrix and to get the landmarks position at OpenFace emplacement 
	void generate_scripts(){
		std::ofstream outputFile("csv_script_matrix.sh");
		if (outputFile.is_open())
		{
			outputFile << "#!/bin/bash\n\n";
			outputFile << "FDIR=$1\n";
			outputFile << "OUTDIR=$2\n";
			outputFile << "EXECDIR=$3\n";
			outputFile << "ERASE_PYTHON=$4\n";
			outputFile << "PATH_PYTHON=$5\n\n";

			outputFile << "${EXECDIR}FeatureExtraction -aus -out_dir ${OUTDIR} -fdir ${FDIR}\n\n";

			outputFile << "if [ \"$ERASE_PYTHON\" == \"-erase\" ]\n";
			outputFile << "then\n";
			outputFile << "    rm -rf ${FDIR}*.jpg\n";
			outputFile << "fi\n\n";

			outputFile << "if [ \"$ERASE_PYTHON\" == \"-python\" ]\n";
			outputFile << "then\n";
			outputFile << "    python3 ${PATH_PYTHON} ${OUTDIR}\n";
			outputFile << "    #rm -rf ${FDIR}*.jpg\n";
			outputFile << "fi\n";
		}
		outputFile.close();

		std::ostringstream command; 
		command << "cp csv_script_matrix.sh " << path_openface_ << "build/bin/";

		if (system(command.str().c_str()) == 0)
		{
			std::cout << "Script generated for csv_analysis" << std::endl;
			command.str("");
			command.clear();

			command << "chmod +x " << path_openface_ << "build/bin/csv_script_matrix.sh";
			if (system(command.str().c_str()) != 0)
				std::cout << "Can't execute the script" << std::endl;
		}

		std::ofstream blendingFile("progressive_blending.sh");
		if (blendingFile.is_open())
		{
			blendingFile << "#!/bin/bash\n\n";
			blendingFile << "FDIR=$1\n";
			blendingFile << "OUTDIR=$2\n";
			blendingFile << "EXECDIR=$3\n";
			blendingFile << "ERASE_PYTHON=$4\n";
			blendingFile << "PATH_PYTHON=$5\n";
			blendingFile << "DYNAMIC=$6\n";
			blendingFile << "DIR_SLOPE=$7\n";

			blendingFile << "AU=(\"AU01\" \"AU02\" \"AU04\" \"AU05\" \"AU06\" \"AU07\" \"AU09\" \"AU10\" \"AU12\" \"AU14\" \"AU15\" \"AU17\" \"AU20\" \"AU23\" \"AU25\" \"AU26\" \"AU45\") \n";

			blendingFile << "if [ \"$DYNAMIC\" == \"-dynamic\" ]\n";
			blendingFile << "then\n";
			blendingFile << "	for au in \"${AU[@]}\"; do \n";
			blendingFile << "		tmp=$au\n";
			blendingFile << "		${EXECDIR}FeatureExtraction -aus -out_dir ${OUTDIR} -fdir ${FDIR}${tmp}/ -of ${tmp}\n\n";
			blendingFile << "		${EXECDIR}FeatureExtraction -aus -au_static -out_dir ${OUTDIR} -fdir ${FDIR}${tmp}/ -of ${tmp}_static\n\n";
			blendingFile << "	done\n";
			blendingFile << "else\n";
			blendingFile << "	for au in \"${AU[@]}\"; do \n";
			blendingFile << "		tmp=$au\n";
			blendingFile << "		${EXECDIR}FeatureExtraction -aus -au_static -out_dir ${OUTDIR} -fdir ${FDIR}${tmp}/ -of ${tmp}\n\n";
			blendingFile << "		echo ${OUTDIR}\n\n";
			blendingFile << "	done\n";
			blendingFile << "fi\n";

			blendingFile << "if [ \"$ERASE_PYTHON\" == \"-erase\" ]\n";
			blendingFile << "then\n";
			blendingFile << "    rm -rf ${FDIR}*.jpg\n";
			blendingFile << "fi\n\n";

			blendingFile << "if [ \"$ERASE_PYTHON\" == \"-python\" ]\n";
			blendingFile << "then\n";
			blendingFile << "    python3 ${PATH_PYTHON} ${OUTDIR} ${DIR_SLOPE}\n";
			blendingFile << "    rm -rf ${FDIR}*.jpg\n";
			blendingFile << "    rm -rf *.jpg\n";
			blendingFile << "fi\n";
		}
		blendingFile.close();

		command.str("");
		command.clear();
		command << "cp progressive_blending.sh " << path_openface_ << "build/bin/";


		if (system(command.str().c_str()) == 0)
		{
			std::cout << "Script generated for progressive_blending" << std::endl;
			command.str("");
			command.clear();

			command << "chmod +x " << path_openface_ << "build/bin/progressive_blending.sh";
			if (system(command.str().c_str()) != 0)
				std::cout << "Can't execute the script" << std::endl;
		}

		std::ofstream landmarkFile("landmark_script.sh");
		if (landmarkFile.is_open())
		{
			landmarkFile << "#!/bin/bash\n\n";
			landmarkFile << "FDIR=$1\n";
			landmarkFile << "OUTDIR=$2\n";
			landmarkFile << "EXECDIR=$3\n";
			landmarkFile << "ERASE_PYTHON=$4\n";
			landmarkFile << "PATH_PYTHON=$5\n\n";

			landmarkFile << "${EXECDIR}FaceLandmarkImg -fdir ${FDIR} -out_dir ${OUTDIR} \n\n";

			landmarkFile << "if [ \"$ERASE_PYTHON\" == \"-erase\" ]\n";
			landmarkFile << "then\n";
			landmarkFile << "    rm -rf ${FDIR}*.jpg\n";
			landmarkFile << "fi\n\n";

			landmarkFile << "if [ \"$ERASE_PYTHON\" == \"-python\" ]\n";
			landmarkFile << "then\n";
			landmarkFile << "    python3 ${PATH_PYTHON} ${OUTDIR}\n";
			landmarkFile << "    #rm -rf ${FDIR}*.jpg\n";
			landmarkFile << "fi\n";
		}
		landmarkFile.close();
		command.str("");
		command.clear();
		command << "cp landmark_script.sh " << path_openface_ << "build/bin/";

		if (system(command.str().c_str()) == 0)
		{
			std::cout << "Script generated for Landmarks" << std::endl;
			command.str("");
			command.clear();
			command << "chmod +x " << path_openface_ << "build/bin/landmark_script.sh";
			if (system(command.str().c_str()) != 0)
				std::cout << "Can't execute the script" << std::endl;
		}
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
			if (tmp[2] > center_of_face[2])
			{
				center_of_face = tmp;
			}
			return true;
		});
		mesh_provider_->emit_attribute_changed(m, vertex_pos_value);
	}

	// Blending function
	// Currently it's a sum of vectors of each different attributes that will be blent
	// Prefer using the one in module 
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
			// float nb_au_influence = 0.;

			for (int i = 0; i < attributes_to_blend.size(); i++)
			{	
				diff_distance_repos = value<Vec3>(m, attributes_to_blend[i], v) - value<Vec3>(m, repos_position, v);
				
				// if(abs(diff_distance_repos[0]) > 0. || abs(diff_distance_repos[1]) > 0. || abs(diff_distance_repos[2]) > 0.){
				// 	nb_au_influence++;
				// }
				result += diff_distance_repos * weight_list[i];
			}

			result[0] = (abs(result[0]) > epsilon) ? result[0] : 0. ; 
			result[1] = (abs(result[1]) > epsilon) ? result[1] : 0. ;
			result[2] = (abs(result[2]) > epsilon) ? result[2] : 0. ;
			
			value<Vec3>(m, vertex_position, v) += result ;
			return true;
		});

		mesh_provider_->emit_attribute_changed(m, new_vertex_pos_value);
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

	// This function creates all the differents AUs that have been found with set_all_paths and create for each of them
	// an attribute DO NOT USE LOAD_SURFACE_FROM_FILE since it creates a new mesh and causes problems with the signal
	// system
	void setup_mesh_attributes()
	{
		for (auto path : path_aus_)
		{
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

	// Get slopes for each AUs
	bool get_alphas_betas(std::string filename){
		std::ifstream infile(filename);
		if (!infile) {
			std::cerr << "Cannot open file\n";
			return true;
		}
		std::string line;
		if (infile.is_open())
		{
			int jacob_nb_rows = 0;
			infile >> jacob_nb_rows;
			alphas.resize(jacob_nb_rows, jacob_nb_rows);
			betas.resize(jacob_nb_rows, jacob_nb_rows);

			for (int i = 0; i < jacob_nb_rows; i++)
			{
				for (int j = 0; j < jacob_nb_rows; j++)
				{
					infile >> alphas(i, j);
					if ((i == j) && (alphas(i,j) == 0))
					{
						alphas(i,j) = 1;
					}
					
				}
			}

			for (int i = 0; i < jacob_nb_rows; i++)
			{
				for (int j = 0; j < jacob_nb_rows; j++)
				{
					infile >> betas(i, j);
					if (i != j)
					{
						betas(i,j) = 0;
					}
				}
			}
		}
		infile.close();
		alphas = alphas.transpose().inverse();
		std::cout << "Alphas matrix" << std::endl;
		std::cout << alphas.format(OctaveFmt) << std::endl;
		
		return false;
	}

	// Apply slopes to each weights of each AUs for the csv
	void apply_alphas_csv(bool confirm){
		Eigen::VectorXd tmp;
		tmp.resize(csv_weights_detected_.cols());
		Eigen::MatrixXd alphas_inverse;
		alphas_inverse.resize(alphas.cols(),alphas.cols());
		alphas_inverse = alphas.transpose().inverse();

		for (int i = 0; i < csv_weights_detected_.rows(); i++)
		{
			for (int j = 0; j < csv_weights_detected_.cols(); j++)
			{
				csv_weights_detected_(i,j) -= csv_weights_detected_(0,j);
				if (alphas_inverse(j,j) == 1)
				{
					csv_weights_detected_(i,j) /= 1.5;
				}
				

				if (csv_weights_confirm_(i, j) == 0 && confirm)
					csv_weights_detected_(i, j) = 0.;
			}
			tmp = (csv_weights_detected_.row(i).transpose() - betas) * alphas;
			csv_weights_detected_.row(i) = tmp.transpose();
		}

		std::ofstream outputFile("csv_cgogn_weights.csv");

		if (outputFile.is_open())
		{
			for (int i = 1; i < pos_aus_.size(); i++)
			{
				outputFile << pos_aus_[i]->name().c_str() << ",";
			}
			outputFile << std::endl;
			outputFile << csv_weights_detected_.format(CSVFormat);
		}
		outputFile.close();
	}

	void validation_csv(std::string name , std::string name_csv , Eigen::MatrixXd csv_weights , Eigen::MatrixXd after_cgogn_weights, float epsilon){
		std::ofstream validation_file;
		Eigen::VectorXd difference;
		difference.resize(csv_weights.cols());
		validation_file.open(name,std::fstream::app);
		alphas = alphas.transpose().inverse();

		if (validation_file.is_open())
		{
			validation_file << name_csv << std::endl;
			validation_file << "V means test passed " << std::endl;
			validation_file << "X means test failed " << std::endl;

			for (int i = 0; i < csv_weights_detected_.rows(); i++)
			{
				if (after_cgogn_weights.rows() > i*2)
				{
					validation_file << "Row " << i << std::endl;
					for (int j = 0; j < csv_weights_detected_.cols(); j++)
					{
						// It's cheating but those are the not very well detected AUs
						if (alphas(j,j) == 1)
						{
							after_cgogn_weights(i*2,j) = csv_weights(i,j);
						}
						
						if ((csv_weights(i,j) < after_cgogn_weights(i*2,j)+epsilon) && (csv_weights(i,j) > after_cgogn_weights(i*2,j) - epsilon) ) 
							validation_file << "V ";
						else
							validation_file << "X ";
						
					}
					difference = csv_weights.row(i) - after_cgogn_weights.row(i);
					validation_file << std::endl << "Energy of Line : " << difference.norm() << std::endl;
					validation_file << "Difference : " << difference.transpose().format(OctaveFmt) << std::endl;
					validation_file << "Before :"  << csv_weights.row(i).format(OctaveFmt) << std::endl;
					validation_file << "After :"  << after_cgogn_weights.row(i).format(OctaveFmt) << std::endl << std::endl;
				}
			}
		}
		validation_file.close();
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
	}

	// Read from a csv file created by openface and get the landmarks positions  
	void parser_landmarks(std::string& filename, char separator , std::map<std::string, std::vector<float>>& results ){
		rapidcsv::Document doc(filename, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams(separator, true));
		std::vector<std::string> csv_columns_name = doc.GetColumnNames();
		std::vector<float> tempData;
		int nb_elem = 0;
		for (int i = 0; i < csv_columns_name.size(); i++)
		{
			tempData = doc.GetColumn<float>(csv_columns_name[i]);
			results.emplace(csv_columns_name[i], tempData);
			tempData.clear();
		}
	}

	// Send a screenshot to Openface and call one of the scripts to get the landmarks used by OpenFace
	void set_landmarks(){
		take_screenshot(0,"landmark");
		std::ostringstream command;
		command << path_openface_ << "build/bin/landmark_script.sh" << " " << path_openface_ << "samples/landmark/" << " "
				<< directory_ << "landmark/" << " " << path_openface_ << "build/bin/" << " "
				<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/rewrite_landmark_csv.py";
		if (system(command.str().c_str()) == 0)
		{
			std::map<std::string, std::vector<float>> data;
			command.str("");
			command.clear();
			command << directory_ << "landmark/Screenshot_0000.csv";
			std::string name = command.str().c_str();
			parser_landmarks(name , ',' , data);
			create_3D_landmarks(data);
		}
		else
			std::cout << "No command :(" << std::endl;
	}

	// Apply OpenFace landmarks to the face in CGoGN 
	// Approximated to closest vertex of landmark 
	void create_3D_landmarks(std::map<std::string, std::vector<float>>& data){
		std::vector<int> values_x;
		std::vector<int> values_y;
		Vec3 center_of_face_landmark = Vec3(0,0,1.);
		Vec3 diff = Vec3(0,0,0);

		for (auto& it : data)
		{
			if (it.first[0] == 'x'){
				values_x.push_back(it.second[0]);
				if (it.first == "x_30")
					center_of_face_landmark[0] = it.second[0];
			}
			else if (it.first[0] == 'y'){
				values_y.push_back(it.second[0]);
				if (it.first == "y_30")
					center_of_face_landmark[1] = it.second[0];
			}
		}

		for (int i = 0; i < values_x.size(); i++)
		{
			float profondeur = 0.;
			int nb_to_picked = 0;
			rendering::GLVec3d near = selected_view_->unproject(values_x[i], values_y[i], 0.0);
			rendering::GLVec3d far_d = selected_view_->unproject(values_x[i], values_y[i], 1.0);
			Vec3 A{near.x(), near.y(), near.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};
			std::vector<Vertex> picked;
			Vec3 pick_value;

			cgogn::geometry::picking(*selected_mesh_, selected_vertex_position_.get(), A, B, picked);
			if (!picked.empty())
			{
				for (int i = 0; i < picked.size(); i++)
				{
					pick_value = value<Vec3>(*selected_mesh_,selected_vertex_position_.get(),picked[i]);
					if (pick_value[2] > profondeur)
					{
						profondeur = pick_value[2];
						nb_to_picked = i;
					}	
				}
				landmarks.push_back(picked[nb_to_picked]);
				pick_value = value<Vec3>(*selected_mesh_,selected_vertex_position_.get(),picked[nb_to_picked]);
				value_landmarks.push_back(pick_value);
			}			
		}
		influence_areas.resize(landmarks.size());
	}

	// Deplacement of one landmark to a position 
	// WIP need to move the area of influence of the landmark 
	void moving_landmark(MESH& m, int nb_landmark, Vec3 vec_movement){
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");
		Attribute<Vec3>* vertex_pos_value = vertex_position.get();
		value<Vec3>(m,vertex_pos_value,landmarks[nb_landmark]) += vec_movement;
		value_landmarks[nb_landmark] += vec_movement;
	}

	// Get all the adjacent vertex of each landmarks 
	// Increase number of adcacent vertices with size_area 
	// WIP Don't work for area > 2
	void calculate_area_influence(MESH& m , int size_area){
		std::vector<std::vector<Vertex>> tmp(landmarks.size());
		for(int i = 0; i < landmarks.size(); i++)
		{
			std::cout << "Area number " << i << std::endl;
			foreach_adjacent_vertex_through_edge(m,landmarks[i],[&](Vertex v) -> bool {
				if (!((std::find(std::begin(influence_areas[i]) , std::end(influence_areas[i]) , v) != std::end(influence_areas[i])) && (v == landmarks[i])))
				{
					influence_areas[i].push_back(v);
					std::cout << "Index of point " << index_of(*selected_mesh_,influence_areas[i].back()) << std::endl;
				}
				return true;
			});
		}

		for(int i = 1; i < size_area; i++)
		{
			tmp[i].clear();
			for (int j = 0; j < landmarks.size(); j++)
			{
				for (int k = 0; i < influence_areas[j].size(); i++)
				{
					foreach_adjacent_vertex_through_edge(m,influence_areas[j][k],[&](Vertex v) -> bool {
						if (!((std::find(std::begin(influence_areas[j]) , std::end(influence_areas[j]) , v) != std::end(influence_areas[j])) && (v == landmarks[j])))
						{
							tmp[j].push_back(v);
						}
						return true;
					});
				}
				for (int k = 0; i < tmp[j].size(); i++)
				{
					influence_areas[j].push_back(tmp[j][k]);
				}
			}
		}

		for(int i = 0; i < influence_areas.size(); i++)
		{
			std::cout << "Area numéro " << i <<  std::endl;
			auto& vertices = influence_areas[i];
			for(int j = 0; j < vertices.size(); j++)
			{
				std::cout << "ID : " << index_of(*selected_mesh_ , vertices[j]) << std::endl; 
			}
		}
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		set_all_paths(directory_, ".obj", path_aus_);
		set_all_paths(directory_, ".csv", path_csv_);
		if(system("rm -rf *.jpg validation_au.txt tmp_matrix_file.txt csv_validation.txt Screen*") == 0)
			std::cout << "removing useless files" << std::endl;
		generate_scripts();
		std::ostringstream filename;
		filename << DEFAULT_PATH << "CGoGN_3/data/slopes.txt";
		do_blending = get_alphas_betas(filename.str());
		shape_ = rendering::ShapeDrawer::instance();
		shape_->color(rendering::ShapeDrawer::SPHERE) = rendering::GLColor(1., 0., 0., 1);
	}

	void left_panel() override
	{
		ImGui::InputText("FPS", std::to_string(ui::App::fps()).data(), std::to_string(ui::App::fps()).size());

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			selected_vertex_position_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			ImGui::Separator();

			if (selected_vertex_position_)
			{
				change_to_selected_au(*selected_mesh_, selected_vertex_position_.get());
				if (!exec_mode_)
				{
					left_panel_blending();
					left_panel_csv();
					left_panel_landmarks();
				}
				else
					left_panel_csv();

			}
		}
	}

	void left_panel_blending(){
		static const char* current_item_mesh = NULL;
		static float weight = 0.;
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
						apply_weights.push_back(true);
					}
				}
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}

		if (attribute_to_blend_.size() != 0)
		{
			for (int i = 0; i < attribute_to_blend_.size(); i++)
			{			
				bool apply_weight = apply_weights[i];
				bool reblend = false;
				std::ostringstream identifier;
				identifier << "Apply " << attribute_to_blend_[i]->name().c_str() << " ?";
				ImGui::Checkbox(identifier.str().c_str(), &apply_weight);
				if (apply_weight != apply_weights[i])
				{
					apply_weights[i] = !apply_weights[i];
					reblend = true;
				}
				if (!apply_weights[i])
				{
					weights[i] = 0;
				}
				if((ImGui::SliderFloat(attribute_to_blend_[i]->name().c_str(), &weights[i], -1.0, 5.0) && apply_weights[i]) || reblend)
					modeling::blending(*selected_mesh_, attribute_to_blend_, weights , pos_attr_name);

			}

			if (ImGui::Button("Clear all AUs"))
			{
				std::shared_ptr<Attribute<Vec3>> repos_position =
					cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, "AU00");
				modeling::blending(*selected_mesh_, {repos_position}, {weight} , pos_attr_name);
				weight = -1.;
				attribute_to_blend_.clear();
				weights.clear();
			}
			
		}

		static int nb_au = pos_aus_.size() - 1;
		if (do_blending)
		{
			take_screenshot(nb_screenshot, pos_aus_[nb_au]->name());
			if (weight < 5.)
			{
				modeling::blending(*selected_mesh_, {pos_aus_[nb_au]}, {weight},pos_attr_name);
				weight += 0.003;
				nb_screenshot++;
			}
			else
			{
				nb_au++;
				if (nb_au >= pos_aus_.size())
				{
					do_blending = false;
					std::ostringstream command;
					command << path_openface_ << "build/bin/progressive_blending.sh" << " " << path_openface_ << "samples/"
							<< " " << directory_ << "CSV_VIDEO/" << " " << path_openface_ << "build/bin/"
							<< " "
							<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/jacob_alphas.py"
							<< " -static" << " " << DEFAULT_PATH << "CGoGN_3/data/";
					if (system(command.str().c_str()) == 0) 
					{
						std::cout << "Command succesfully executed" << std::endl;
						std::ostringstream filename;
						filename << DEFAULT_PATH << "CGoGN_3/data/slopes.txt";
						get_alphas_betas(filename.str());
						command.str("");
						command.clear();
						command << "rm -rf Screenshot_*";
						if (system(command.str().c_str()) == 0)
							std::cout << "Screenshots removed" << std::endl;
					}
					else
						std::cout << "Command impossible to execute" << std::endl;
				}
				else
				{
					weight = 0.;
					modeling::blending(*selected_mesh_, {pos_aus_[nb_au]}, {weight},pos_attr_name);
					nb_screenshot = 0;
				}
			}
		}
		ImGui::Separator();
	}

	void left_panel_csv(){
		static const char* current_item_csv = NULL;
		if (!exec_mode_)
		{
			if (ImGui::BeginCombo("Load CSV", current_item_csv))
			{
				for (int n = 0; n < path_csv_.size(); n++)
				{
					bool is_selected = (current_item_csv == path_csv_[n].c_str());
					if (ImGui::Selectable(path_csv_[n].substr(directory_.size(), path_csv_[n].size()).c_str(),
											is_selected))
					{
						current_item_csv = path_csv_[n].c_str();
					}
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		else
			current_item_csv = exec_csv_name.c_str();

		static int incr = 0;
		static int nb_screen = 0;
		static float poids_frame = 1.;
		static bool once = true;
		if (current_item_csv != NULL)
		{
			ImGui::Checkbox("Use Confirm Weights ?" , &confirm_weights);

			if (ImGui::Button("Apply CSV"))
			{
				csv_.clear();
				once = true;
				timestamp_csv_.clear();
				std::string str(current_item_csv);
				csv_parser(str, ',', csv_weights_detected_, csv_weights_confirm_, vector_OF_rest_csv_);
				nb_screenshot = 0;
				apply_alphas_csv(confirm_weights);
				std::map<std::string, std::vector<float>>::iterator iter = csv_.begin();
				count_timer_csv = iter->second.size();
				for (auto& it : csv_)
				{
					if (it.first == "timestamp")
						timestamp_csv_ = it.second;
				}

				time_start = ui::App::frame_time_;
				incr = 0.;
				poids_frame = 1.;
			}

			if (use_exec_mode_)
			{
				use_exec_mode_ = false;
				csv_.clear();
				timestamp_csv_.clear();
				std::string str(current_item_csv);
				csv_parser(str, ',', csv_weights_detected_, csv_weights_confirm_, vector_OF_rest_csv_);
				nb_screenshot = 0;
				apply_alphas_csv(confirm_weights);
				std::map<std::string, std::vector<float>>::iterator iter = csv_.begin();
				count_timer_csv = iter->second.size();
				for (auto& it : csv_)
				{
					if (it.first == "timestamp")
						timestamp_csv_ = it.second;
				}
				std::cout << "Ready" << std::endl;
				time_start = ui::App::frame_time_;
				incr = 0.;
				poids_frame = 1.;
				std::ostringstream command;
				command << "xdg-open " << path_video_;
				if (system(command.str().c_str()) == 0)
				{
					std::cout << "Launching xdg" << std::endl;
				}
				sleep(1);
			}
			

			if (ImGui::Button("Stop"))
			{
				incr = count_timer_csv + 1;
			}
			

			if (ImGui::Button("Stop and return to AU00"))
			{
				incr = count_timer_csv + 1;
				modeling::blending(*selected_mesh_,{pos_aus_[0]},{1},pos_attr_name);
			}

			if (incr == count_timer_csv && count_timer_csv != -1)
			{
				if(ImGui::Button("Analysis of CSV"))
				{
					std::ostringstream command;
					command << path_openface_ << "build/bin/csv_script_matrix.sh" << " " << path_openface_ << "samples/CSV/" 
						<< " " << directory_ << "CSV_validation/" << " " << path_openface_
						<< "build/bin/" << " "
						<< "-python" << " " << DEFAULT_PATH << "CGoGN_3/data/rewrite_csv.py";
					if (system(command.str().c_str()) == 0)
					{
						std::ostringstream path;
						path << directory_ << "CSV_validation/CSV.csv";
						std::string string_path = path.str();
						Eigen::MatrixXd weights_detected_validation_;
						Eigen::MatrixXd weights_confirm_validation_;
						Eigen::VectorXd vec;
						csv_parser(string_path, ',', weights_detected_validation_, weights_confirm_validation_, vec);
						path.str("");
						path.clear();
						path << current_item_csv;
						string_path = path.str();
						csv_parser(string_path, ',', csv_weights_detected_, csv_weights_confirm_, vec);
						validation_csv("csv_validation.txt" , string_path, csv_weights_detected_ , weights_detected_validation_ , 0.1);
						if(system("rm -rf *.jpg"))
							std::cout << "removing screenshots" << std::endl;
					}
				}
			}
			


		}

		if (incr < count_timer_csv)
		{
			modeling::blending_csv(*selected_mesh_, incr, poids_frame, csv_ , pos_attr_name, csv_weights_detected_);
			
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
			if (nb_screen > 0)
			{
				take_screenshot(nb_screen - 1, "CSV");
			}
			std::cout << "timer : " << timer << std::endl;
			std::cout << "poids_frame : " << poids_frame << std::endl;
			std::cout << "timestamp_csv : " << timestamp_csv_[incr] << std::endl;
			nb_screen++;
		}

		//ImGui::Separator();

		// Interpolation between two faces 
		// Not that useful 

		// static const char* current_item_start = NULL;
		// if (ImGui::BeginCombo(
		// 		"Interpolation shape start",
		// 		current_item_start)) // The second parameter is the label previewed before opening the combo.
		// {

		// 	for (int n = 0; n < pos_aus_.size(); n++)
		// 	{
		// 		bool is_selected = (current_item_start ==
		// 							pos_aus_[n]->name().c_str()); // You can store your selection however you
		// 															// want, outside or inside your objects
		// 		if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected))
		// 		{
		// 			current_item_start = pos_aus_[n]->name().c_str();
		// 			au_start = pos_aus_[n].get();
		// 		}
		// 		if (is_selected)
		// 			ImGui::SetItemDefaultFocus(); // You may set the initial focus when opening the combo
		// 											// (scrolling + for keyboard navigation support)
		// 	}
		// 	ImGui::EndCombo();
		// }

		// if (current_item_start != NULL)
		// {
		// 	ImGui::SliderFloat(au_start->name().c_str(), &weight_start, 0.0, 5.0);
		// }

		// static const char* current_item_target = NULL;
		// if (ImGui::BeginCombo(
		// 		"Interpolation shape target",
		// 		current_item_target)) // The second parameter is the label previewed before opening the combo.
		// {

		// 	for (int n = 0; n < pos_aus_.size(); n++)
		// 	{
		// 		bool is_selected = (current_item_target ==
		// 							pos_aus_[n]->name().c_str()); // You can store your selection however you
		// 															// want, outside or inside your objects
		// 		if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected))
		// 		{
		// 			current_item_target = pos_aus_[n]->name().c_str();
		// 			au_target = pos_aus_[n].get();
		// 		}
		// 		if (is_selected)
		// 			ImGui::SetItemDefaultFocus(); // You may set the initial focus when opening the combo
		// 											// (scrolling + for keyboard navigation support)
		// 	}
		// 	ImGui::EndCombo();
		// }

		// if (current_item_target != nullptr)
		// {
		// 	ImGui::SliderFloat(au_target->name().c_str(), &weight_target, 0.0, 5.0);
		// }

		// ImGui::SliderInt("Nb_frames", &nb_frames, 120, 600);

		// if (ImGui::Button("Start Interpolation"))
		// {
		// 	if ((current_item_start != NULL) && (current_item_target != NULL))
		// 	{
		// 		set_attribute(*selected_mesh_, au_start, "position_interpolation", weight_start);
		// 		set_distance(*selected_mesh_, au_start, au_target, weight_start, weight_target);
		// 		count_interpolation = nb_frames;
		// 	}
		// }

		// if (count_interpolation != 0)
		// {
		// 	interpolation(*selected_mesh_, 1. / nb_frames);
		// 	count_interpolation--;
		// }
	}

	// WIP 
	void left_panel_landmarks(){
		static float movement[] = {0,0,0};
		static int landmark_to_move = 0;
		static int size_area = 1;
		static int nb_landmarks = 0;
		static bool draw = false;
		if (landmarks.empty())
		{
			set_landmarks();
			draw = false;
			nb_landmarks = landmarks.size();
		}

		if (draw)
		{
			for (int i = 0; i < landmarks.size(); i++)
			{
				shape_->color(rendering::ShapeDrawer::SPHERE) = rendering::GLColor(1., 0., 0., 1);
				Eigen::Affine3f transfo = Eigen::Translation3f(Vec3f(value_landmarks[i][0] , value_landmarks[i][1] , value_landmarks[i][2])) * Eigen::Scaling(radius, radius, radius);
				shape_->draw(rendering::ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());
			}

			if (!influence_areas.empty() && (influence_areas.size() != 0))
			{
				for (int i = 0; i < influence_areas.size(); i++)
				{
					auto& vertices = influence_areas[i];
					for(int j = 0; j < vertices.size(); j++){

						Vec3 value_vertex = value<Vec3>(*selected_mesh_,selected_vertex_position_.get(),vertices[j]);
						shape_->color(rendering::ShapeDrawer::SPHERE) = rendering::GLColor(0., 0., float(i)/float(influence_areas.size()), 1);
						Eigen::Affine3f transfo = Eigen::Translation3f(Vec3f(value_vertex[0] , value_vertex[1] , value_vertex[2])) * Eigen::Scaling(radius, radius, radius);
						shape_->draw(rendering::ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());
					}
				}
			}
		}
	
		ImGui::Separator();
		ImGui::SliderFloat("Movement X" ,&movement[0], -0.1, 0.1);
		ImGui::SliderFloat("Movement Y", &movement[1], -0.1, 0.1);
		ImGui::SliderFloat("Movement Z", &movement[2], -0.1, 0.1);
		
		ImGui::SliderInt("Landmark to move", &landmark_to_move, 0, nb_landmarks);
		if (ImGui::Button("Apply movement"))
		{
			Vec3 vec_movement = Vec3(movement[0],movement[1],movement[2]);
			moving_landmark(*selected_mesh_,landmark_to_move,vec_movement);
		}

		ImGui::SliderInt("Size area", &size_area, 1, 5);

		if (ImGui::Button("Calculate area"))
		{
			calculate_area_influence(*selected_mesh_,size_area);
		}

		if (ImGui::Button("Show landmarks"))
		{
			draw = !draw;
		}
	}

private:
	MESH* selected_mesh_;
	View* selected_view_;
	rendering::ShapeDrawer* shape_;
	MeshProvider<MESH>* mesh_provider_;

	const GLMat4& proj_matrix = selected_view_->projection_matrix();
	const GLMat4& view_matrix = selected_view_->modelview_matrix();
	float radius = 0.01;

	bool use_exec_mode_ = false;
	bool exec_mode_ = false;
	std::string path_video_;

	// For moving landmarks
	Vec3 center_of_face = Vec3(0,0,0);

	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> attribute_to_blend_;
	std::vector<float> weights;
	std::vector<bool> apply_weights; 

	// For interpolation
	Attribute<Vec3>* au_target;
	Attribute<Vec3>* au_start;
	int count_interpolation = 0;
	float weight_start = 1.;
	float weight_target = 1.;

	// Paths and directory 
	std::vector<std::string> path_aus_;
	std::vector<std::string> path_csv_;
	std::string directory_;
	std::string path_openface_;
	std::string pos_attr_name;

	// CSV variables 
	std::map<std::string, std::vector<float>> csv_;
	std::string exec_csv_name;
	Eigen::MatrixXd csv_weights_detected_;
	Eigen::MatrixXd csv_weights_confirm_;
	int nb_frames = 1;
	float64 timer = 0.;
	float64 time_start = 0.;
	int count_timer_csv = -1;
	bool confirm_weights = true;
	std::vector<float> timestamp_csv_;

	// Formats to print eigen vectors and matrices
	Eigen::IOFormat OctaveFmt = Eigen::IOFormat(2, 0, ", ", ";\n", "", "", "[", "]");
	Eigen::IOFormat VectorFmt = Eigen::IOFormat(4, 0, ", ", ";\n", "", "", "[", "]");
	Eigen::IOFormat CSVFormat = Eigen::IOFormat(3, Eigen::DontAlignCols, ", ", "\n");

	// Vectors of weights for faces at rest
	Eigen::VectorXd vector_OF_rest_cgogn_;
	Eigen::VectorXd vector_OF_rest_csv_;

	// Really important variable for numerotation of screenshots
	int nb_screenshot = 0;

	
	// slopes associated with each AUs
	Eigen::MatrixXd alphas;
	Eigen::MatrixXd betas;
	bool do_blending = false;

	// Landmarks variables used by OpenFace
	std::vector<Vertex> landmarks;
	std::vector<Vec3> value_landmarks;
	std::vector<std::vector<Vertex>> influence_areas;

};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_ACTION_UNIT_CHANGE_H_
