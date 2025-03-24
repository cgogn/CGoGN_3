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
#include <sstream>
#include <string>

#define PATH_OF "../../../../OpenFace/OpenFace/build/bin/csv_script_matrix.sh"
#define DEFAULT_PATH CGOGN_STR(CGOGN_DATA_PATH) "../../"

namespace fs = std::filesystem;

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;

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
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, vertex_pos_value, v) = value<Vec3>(m, au_position, v);
			return true;
		});
		mesh_provider_->emit_attribute_changed(m, vertex_pos_value);
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

	void setup_csv_matrix()
	{
		int i = 0;
		std::ostringstream command;
		command << PATH_OF << " " << DEFAULT_PATH << "OpenFace/OpenFace/samples/Matrix/" << " " << directory_ << "CSV/" << " " << DEFAULT_PATH << "OpenFace/OpenFace/build/bin/";
		for (auto au : pos_aus_)
		{
			highlight_difference(*selected_mesh_, au);
			blending(*selected_mesh_, au, epsilon);

			std::ostringstream name;
			name << "/home/roger/Documents/stage/CGoGn/CGoGN_3/build/stage/bin/";
			name << "Screenshot_";
			for (int j = 0; j < 4 - std::to_string(i).size(); j++)
			{
				name << "0";
			}
			name << i;
			name << ".jpg";

			std::ostringstream dirname;
			dirname << "/home/roger/Documents/stage/CGoGn/OpenFace/OpenFace/samples/Matrix/" ;

			selected_view_->save_screenshot_name(name.str());
			fs::path sourceFile = name.str().c_str();
			fs::path targetParent = dirname.str().c_str();
			if (!fs::is_directory(targetParent) || !fs::exists(targetParent))
				fs::create_directory(targetParent);
			
			fs::copy(sourceFile, targetParent, fs::copy_options::overwrite_existing);

			if (strcmp(au->name().c_str(), "AU00") == 0)
			{
				if(system(command.str().c_str()) == 0)
				{
					std::ostringstream matrix_csv_path;
					matrix_csv_path << directory_ << "CSV/Matrix.csv";
					std::string path = matrix_csv_path.str();
					csv_parser(path, ',');
					std::vector<int> weights_confirm;
					std::vector<float> weights_detect;
					for (auto& it : csv_)
					{
						if (ends_with(it.first, "_r"))
							weights_detect.push_back(it.second[0]);
						if (ends_with(it.first, "_c"))
							weights_confirm.push_back(it.second[0]);
					}
					for (int i = 0; i < weights_confirm.size(); i++)
					{
						// if (weights_confirm[i] == 0)
						vector_OF_rest_cgogn.push_back(weights_detect[i]);
						// else
						// 	vector_OF_rest_cgogn.push_back(0);

						std::cout << vector_OF_rest_cgogn[i] << " ";
					}
					std::cout << std::endl;
				}
				else
					std::cout << "Command invalid" << std::endl;
			}
			else
				i++;
		}
		
		if(system(command.str().c_str()) == 0)
		{
			if(system("rm -rf *.jpg") == 0)
				std::cout << "Erasing screenshots" << std::endl;
			else
				std::cout << "Command invalid" << std::endl;

			std::ostringstream matrix_csv_path;
			matrix_csv_path << directory_ << "CSV/Matrix.csv";
			// path_csv_.push_back(matrix_csv_path.str());
			std::string path = matrix_csv_path.str();
			csv_parser(path, ',');
			get_matrix_csv(path);
		}
		else
			std::cout << "Command invalid" << std::endl;
	}

	void get_matrix_csv(std::string path_csv)
	{
		std::vector<std::vector<float>> weights_confirm;
		std::vector<std::vector<float>> weights_detect;
		int incr = 0;
		for (auto& it : csv_)
		{	
			std::cout << it.first.c_str() << std::endl;
			if (ends_with(it.first, "_r"))
				weights_detect.push_back(it.second);
			if (ends_with(it.first, "_c"))
				weights_confirm.push_back(it.second);
		}
		std::cout << weights_detect.size() << weights_detect[0].size() << std::endl;
		for (int i = 0; i < weights_confirm.size(); i++)
		{
			std::vector<float> tmp;
			for (int j = 0; j < weights_confirm[i].size(); j++)
			{
				
				if (weights_confirm[i][j] == 0.){
					weights_detect[i][j] = 0.;
				}
				tmp.push_back(0.);
			}
			matrix_jacob.push_back(tmp);
		}

		for (int i = 0; i < weights_detect[0].size(); i++)
		{
			std::cout << "Line number : " << i << std::endl;
			for (int j = 0; j < weights_detect.size(); j++)
			{
				float element_jacob = (((vector_OF_rest_cgogn[j] + (weights_detect[j][i] * epsilon)) - vector_OF_rest_cgogn[j] ) / epsilon);
				std::cout << ((element_jacob > epsilon) ? element_jacob : 0) << " "; 
				matrix_jacob[j][i] = ((element_jacob > epsilon) ? element_jacob : 0);
				std::cout << matrix_jacob[j][i] << " ";
			}
			std::cout << std::endl;
		}
	}

	// Parse a csv using a filename and the separator of the csv
	void csv_parser(std::string& filename, char separator)
	{
		rapidcsv::Document doc(filename, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams(separator, true));
		std::ofstream outputFile("test.txt"); // Open/create a file named "test.txt" for writing
		std::vector<std::string> csv_columns_name = doc.GetColumnNames();
		std::vector<float> tempData;
		csv_.clear();
		for (int i = 0; i < csv_columns_name.size(); i++)
		{
			if (strcmp(csv_columns_name[i].c_str() , "AU28_c") != 0)
			{
				std::cout << csv_columns_name[i].c_str() << std::endl;
				tempData = doc.GetColumn<float>(csv_columns_name[i]);
				for (int j = 0; j < tempData.size(); j++)
				{
					if (outputFile.is_open())
					{ // Check if the file was successfully opened
						// Write some text into the file
						outputFile << tempData[i];
						outputFile << ";";
						// Close the file
					}
					else
					{
						std::cout << "Failed to create the file."
								<< std::endl; // Display an error message if file creation failed
					}
				}
				outputFile << tempData.size();
				outputFile << "\n";
				csv_.emplace(csv_columns_name[i], tempData);
				tempData.clear();
			}
		}
		std::cout << "Text has been written to the file." << std::endl; // Display a success message
		outputFile.close();												// Close the file after writing

		std::ofstream outputFile2("test2.txt"); // Open/create a file named "test2.txt" for writing
		int i = 0 ;
		for (auto& it : csv_)
		{
			std::cout << i << " ";
			i++;
			outputFile2 << it.first;
			outputFile2 << ";";
			for (auto& i : it.second)
			{
				if (outputFile2.is_open())
				{ // Check if the file was successfully opened
					// Write some text into the file
					outputFile2 << i;
					outputFile2 << ";";
					// Close the file
				}
				else
					std::cout << "Failed to create the file."
							  << std::endl; // Display an error message if file creation failed
			}

			outputFile2 << it.second.size();
			outputFile2 << "\n";
		}
		std::cout << std::endl;
		std::cout << "Text has been written to the file." << std::endl; // Display a success message
		outputFile2.close();											// Close the file after writing
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
	void blending(MESH& m, std::shared_ptr<Attribute<Vec3>> attribute_to_blend, float weight)
	{
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");
		std::shared_ptr<Attribute<Vec3>> color = cgogn::get_attribute<Vec3, Vertex>(m, "color");
		std::shared_ptr<Attribute<Vec3>> repos_position = cgogn::get_attribute<Vec3, Vertex>(m, "AU00");
		Attribute<Vec3>* new_vertex_pos_value = vertex_position.get();
		Vec3 diff_distance_repos = Vec3(0, 0, 0);

		if (attribute_to_blend->name().c_str() == repos_position->name().c_str())
		{
			parallel_foreach_cell(m, [&](Vertex v) -> bool {
				value<Vec3>(m, vertex_position, v) = value<Vec3>(m, repos_position, v);
				return true;
			});
		}
		else
		{
			parallel_foreach_cell(m, [&](Vertex v) -> bool {
				if (value<Vec3>(m, color, v) == GREEN)
				{
					diff_distance_repos = value<Vec3>(m, attribute_to_blend, v) - value<Vec3>(m, repos_position, v);
					diff_distance_repos = diff_distance_repos * weight;
					value<Vec3>(m, vertex_position, v) = value<Vec3>(m, repos_position, v) + diff_distance_repos;
				}
				return true;
			});
		}
		mesh_provider_->emit_attribute_changed(m, new_vertex_pos_value);
	}

	// Function called each frame after clicking on Apply CSV
	// Setup the differents AUs and weights to calculate the frames
	void blending_csv(std::vector<float> weights, int incr, float poids_frame)
	{
		std::vector<int> aus_confirm;
		std::vector<int> aus_confirm_next_frame;
		std::vector<float> weights_frame;
		std::vector<float> weights_next_frame;
		std::vector<std::shared_ptr<Attribute<Vec3>>> attributes_csv;

		for (auto& it : csv_)
		{
			if (ends_with(it.first, "_r"))
			{
				attributes_csv.push_back(
					cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, it.first.substr(0, it.first.size() - 2)));
				if (incr + 1 < count_timer_csv)
				{
					weights_frame.push_back(it.second[incr]);
					weights_next_frame.push_back(it.second[incr + 1]);
				}
				else
				{
					weights_frame.push_back(it.second[incr - 1]);
					weights_next_frame.push_back(it.second[incr]);
				}
				if (it.second[incr] > 0.01)
					std::cout << it.first << " poids : " << it.second[incr] << std::endl;
			}
			if (ends_with(it.first, "_c"))
			{
				if (incr + 1 < count_timer_csv)
				{
					aus_confirm.push_back(it.second[incr]);
					aus_confirm_next_frame.push_back(it.second[incr + 1]);
				}
				else
				{
					aus_confirm.push_back(it.second[incr - 1]);
					aus_confirm_next_frame.push_back(it.second[incr]);
				}
			}
		}
		for (int i = 0; i < attributes_csv.size(); i++)
		{
			if (aus_confirm[i] == 0)
				weights_frame[i] = 0.;
			if (aus_confirm_next_frame[i] == 0)
				weights_next_frame[i] = 0.;

			weights.push_back((weights_frame[i] * (1 - poids_frame)) + (weights_next_frame[i] * poids_frame));
			highlight_difference(*selected_mesh_, attributes_csv[i]);
			blending(*selected_mesh_, attributes_csv[i], weights[i]);
		}
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
	}

	void left_panel() override
	{
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
						current_item_mesh)) // The second parameter is the label previewed before opening the combo.
				{

					for (int n = 0; n < pos_aus_.size(); n++)
					{
						bool is_selected = (current_item_mesh ==
											pos_aus_[n]->name().c_str()); // You can store your selection however you
																		  // want, outside or inside your objects
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
							ImGui::SetItemDefaultFocus(); // You may set the initial focus when opening the combo
														  // (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}

				if (attribute_to_blend_.size() != 0)
				{
					for (int i = 0; i < attribute_to_blend_.size(); i++)
					{
						ImGui::SliderFloat(attribute_to_blend_[i]->name().c_str(), &weights[i], -1.0, 2.0);
					}
				}
				static bool start = false;
				if (ImGui::Button("Blend"))
				{
					std::ostringstream new_attribute_name;
					blending(*selected_mesh_, pos_aus_[0], 1.);
					for (int i = 0; i < attribute_to_blend_.size(); i++)
					{
						new_attribute_name << attribute_to_blend_[i]->name().c_str() << 'w' << weights[i] << '+';
						highlight_difference(*selected_mesh_, attribute_to_blend_[i]);
						blending(*selected_mesh_, attribute_to_blend_[i], weights[i]);
					}
					std::shared_ptr<Attribute<Vec3>> new_attribute = add_attribute<Vec3, Vertex>(
						*selected_mesh_,
						new_attribute_name.str().substr(0, new_attribute_name.str().size() - 1).c_str());

					// highlight_difference(*selected_mesh_,attribute_to_blend_[0].get());
					// start = true;
					attribute_to_blend_.clear();
					weights.clear();
				}

				if (ImGui::Button("Blend progressif"))
				{
					highlight_difference(*selected_mesh_, pos_aus_[1]);
					start = true;
					// attribute_to_blend_.clear();
					// weights.clear();
				}

				if (ImGui::Button("Clear"))
				{
					std::shared_ptr<Attribute<Vec3>> repos_position =
						cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, "AU00");
					blending(*selected_mesh_, repos_position, weight);
					start = false;
					weight = -1.;
					attribute_to_blend_.clear();
				}

				if (ImGui::Button("Stop"))
				{
					start = false;
				}

				static int i = 0;
				static int nb_au = 1;
				if (start)
				{
					if (weight < 5.)
					{
						blending(*selected_mesh_, pos_aus_[nb_au], weight);
						weight += 0.003;
					}
					else
					{
						if (nb_au >= pos_aus_.size())
						{
							start = false;
						}
						else
						{
							nb_au++;
							weight = 0.;
							blending(*selected_mesh_, pos_aus_[0], weight);
							highlight_difference(*selected_mesh_, pos_aus_[nb_au]);
							blending(*selected_mesh_, pos_aus_[nb_au], weight);
							i = 0;
						}
					}
				}

				if (ImGui::Button("Screenshot"))
				{
					std::ostringstream name;
					name << "/home/roger/Documents/stage/CGoGn/CGoGN_3/build/stage/bin/";
					name << "Screenshot_" << i << ".png";
					std::cout << name.str().c_str() << std::endl;

					std::ostringstream dirname;
					dirname << "/home/roger/Documents/stage/CGoGn/OpenFace/OpenFace/samples/";
					dirname << selected_vertex_position_->name().c_str() << "/";
					i++;
					std::cout << dirname.str().c_str() << std::endl;

					selected_view_->save_screenshot_name(name.str());
					fs::path sourceFile = name.str().c_str();
					fs::path targetParent = dirname.str().c_str();
					if (!fs::is_directory(targetParent) || !fs::exists(targetParent))
						fs::create_directory(targetParent);
					fs::copy(sourceFile, targetParent, fs::copy_options::overwrite_existing);
				}

				static const char* current_item_csv = NULL;
				if (ImGui::BeginCombo(
						"Load CSV",
						current_item_csv)) // The second parameter is the label previewed before opening the combo.
				{
					for (int n = 0; n < path_csv_.size(); n++)
					{
						bool is_selected =
							(current_item_csv == path_csv_[n].c_str()); // You can store your selection however you
																		// want, outside or inside your objects
						if (ImGui::Selectable(path_csv_[n].substr(directory_.size(), path_csv_[n].size()).c_str(),
											  is_selected))
						{
							current_item_csv = path_csv_[n].c_str();
							csv_.clear();
							timestamp_csv_.clear();
							csv_parser(path_csv_[n], ',');
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus(); // You may set the initial focus when opening the combo
														  // (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}

				static int incr = 0.;
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
					std::cout << "timer : " << timer << std::endl;
					std::cout << "poids_frame : " << poids_frame << std::endl;
					std::cout << "timestamp_csv : " << timestamp_csv_[incr] << std::endl;

					weights.clear();
				}

				ImGui::Separator();

				ImGui::InputText("FPS", std::to_string(ui::App::fps()).data(), std::to_string(ui::App::fps()).size());

				ImGui::InputText("Time since last Frame", std::to_string(ui::App::frame_time_).data(),
								 std::to_string(ui::App::frame_time_).size());

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
					std::ostringstream name;
					name << "/home/roger/Documents/stage/CGoGn/CGoGN_3/build/stage/bin/";

					int nb_num = 4;

					name << "Screenshot_";
					for (int j = 0; j < nb_num - std::to_string(i).size(); j++)
					{
						name << "0";
					}
					name << i;
					name << ".jpg";

					std::cout << name.str().c_str() << std::endl;

					std::ostringstream dirname;
					dirname << "/home/roger/Documents/stage/CGoGn/OpenFace/OpenFace/samples/";
					dirname << pos_aus_[nb_au]->name().c_str() << "/";
					i++;
					std::cout << dirname.str().c_str() << std::endl;

					selected_view_->save_screenshot_name(name.str());
					fs::path sourceFile = name.str().c_str();
					fs::path targetParent = dirname.str().c_str();
					if (!fs::is_directory(targetParent) || !fs::exists(targetParent))
						fs::create_directory(targetParent);
					fs::copy(sourceFile, targetParent, fs::copy_options::overwrite_existing);
				}

				static bool matrix = true;
				if (matrix)
				{
				 	setup_csv_matrix();
					blending(*selected_mesh_,pos_aus_[0],1.);
					matrix = false;
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
	// std::shared_ptr<Attribute<Vec3>> selected_vertex_normal_;

	Attribute<Vec3>* au_target;
	Attribute<Vec3>* au_start;

	std::vector<std::string> path_aus_;
	std::vector<std::string> path_csv_;
	std::string directory_;
	std::map<std::string, std::vector<float>> csv_;

	std::vector<std::vector<float>> matrix_jacob;
	std::vector<float> weights;
	std::vector<float> timestamp_csv_;
	std::vector<float> vector_OF_rest_cgogn;
	std::vector<float> vector_OF_rest_csv;

	int nb_frames = 1;
	float64 timer = 0.;
	float64 time_start = 0.;
	int count_interpolation = 0;
	int count_timer_csv = 0;
	float weight_start = 1.;
	float weight_target = 1.;
	float epsilon = 0.5;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_ACTION_UNIT_CHANGE_H_
