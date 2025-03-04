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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <cstring> 
#include <filesystem>





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
class ActionUnitChange: public ViewModule
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
		  selected_view_(app.current_view()), selected_mesh_(nullptr) , selected_vertex_position_(nullptr)
	{
	}
	~ActionUnitChange()
	{
	}

public:


	static bool ends_with(const std::string& str, const std::string& suffix)
	{
		return str.size() >= suffix.size() && str.compare(str.size()-suffix.size(), suffix.size(), suffix) == 0;
	}

	void set_directory(std::string dirname){
		directory_ = dirname;
	}

	void set_mesh(MESH& m , std::shared_ptr<Attribute<Vec3>> vertex_position){
		selected_mesh_ = &m;
		selected_vertex_position_ = vertex_position;
	}

	void set_view(View& v){
		selected_view_ = &v;
	}

	void set_all(std::string root, std::string ext , std::vector<std::string>& paths)
	{
		root = root.substr(2, root.size() - 3);
		for (auto &p : fs::recursive_directory_iterator(root))
		{
			if (p.path().extension() == ext)
			{
				paths.push_back(p.path().string());
				std::cout << p.path().string() << std::endl;
			}
		}
		std::sort(paths.begin(),paths.end());
	}

	void set_attribute(MESH&m , Attribute<Vec3>* to_set , std::string attribute_name){
		std::shared_ptr<Attribute<Vec3>> attribute_to_change = cgogn::get_or_add_attribute<Vec3, Vertex>(m, attribute_name.c_str());
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m,attribute_to_change,v) = value<Vec3>(m,to_set,v);
			return true;
		});
	}

	std::vector<std::string> get_path(){
		return path_aus_;
	}

	void change_to_selected_au(MESH& m , Attribute<Vec3>* au_position){
		
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");
		Attribute<Vec3>* vertex_pos_value = vertex_position.get();
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m,vertex_pos_value,v) = value<Vec3>(m,au_position,v);
			return true;
		});
		mesh_provider_->emit_attribute_changed(m,vertex_pos_value);
	}

	void setup_mesh_attributes(){

		for (auto path : path_aus_)
		{
			std::cout << path.substr(path.size() - 8 , path.size() - (path.size() - 8) - 4) << std::endl;
			std::ifstream fp(path.c_str(), std::ios::in);
			if (!fp.good())
			{
				std::cerr << "Error opening file " << path.c_str() << std::endl;
				return;
			}
			std::shared_ptr<Attribute<Vec3>> au_pos = cgogn::add_attribute<Vec3, Vertex>(*selected_mesh_, path.substr(path.size() - 8 , path.size() - (path.size() - 8) - 4));
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
					Vec3 temp = Vec3(x,y,z);
					vec_pos.push_back(temp);
				}
			}while (!ss.eof());


			int incr = 0;
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, "position");
			Vec3 point;
			cgogn::foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
				point = value<Vec3>(*selected_mesh_,vertex_position,v);
				std::cout << "Point numéro : " << incr << std::endl;
				std::cout << "Position des points_au : x = " << vec_pos[incr][0] << " y =  " << vec_pos[incr][1] << " z =  " << vec_pos[incr][2]  << std::endl;
				std::cout << "Position des points_repos : x = " << point[0] << " y =  " << point[1] << " z =  " << point[2]  << std::endl;
				cgogn::value<Vec3>(*selected_mesh_,au_pos,v) = vec_pos[incr];
				incr++;
				return true;
			});

			geometry::rescale(*au_pos,1);
			Vec3 point_norm;
			cgogn::foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
				point = value<Vec3>(*selected_mesh_,vertex_position,v);
				point_norm = value<Vec3>(*selected_mesh_,au_pos,v);
				std::cout << "Point numéro : " << incr << std::endl;
				std::cout << "Position des points_au_normalisé : x = " << point_norm[0] << " y =  " << point_norm[1] << " z =  " << point_norm[2]  << std::endl;
				std::cout << "Position des points_repos : x = " << point[0] << " y =  " << point[1] << " z =  " << point[2]  << std::endl;
				return true;
			});
		}
	}

	void csv_parser(std::string& filename,char separator){
		rapidcsv::Document doc(filename,rapidcsv::LabelParams(0,-1),rapidcsv::SeparatorParams(separator,true));
		std::ofstream outputFile("test.txt");  // Open/create a file named "test.txt" for writing
		std::vector<std::string> csv_columns_name = doc.GetColumnNames();
		std::vector<float> tempData;
		for(int i = 0 ; i < csv_columns_name.size() ; i++){
			tempData = doc.GetColumn<float>(csv_columns_name[i]);
			for (int j = 0; j < tempData.size(); j++)
			{
				if (outputFile.is_open()) {  // Check if the file was successfully opened
					// Write some text into the file
					outputFile << tempData[i];
					outputFile << ";";
					// Close the file
				} else {
					std::cout << "Failed to create the file." << std::endl;  // Display an error message if file creation failed
				}
			}
			outputFile << tempData.size();
			outputFile << "\n";
			csv_.emplace(csv_columns_name[i],tempData);
			tempData.clear();
		}
		std::cout << "Text has been written to the file." << std::endl;  // Display a success message
		outputFile.close();  // Close the file after writing

		std::ofstream outputFile2("test2.txt");  // Open/create a file named "test2.txt" for writing

		for (auto &it : csv_)
		{
			outputFile2 << it.first;
			outputFile2 << ";";
			for (auto &i : it.second)
			{
				if (outputFile2.is_open()) {  // Check if the file was successfully opened
					// Write some text into the file
					outputFile2 << i;
					outputFile2 << ";";
					// Close the file
				} else {
					std::cout << "Failed to create the file." << std::endl;  // Display an error message if file creation failed
				}
			}

			outputFile2 << it.second.size();
			outputFile2 << "\n";
		}
		std::cout << "Text has been written to the file." << std::endl;  // Display a success message
		outputFile2.close();  // Close the file after writing
	}

	void set_to_blue(MESH& m){
		std::shared_ptr<Attribute<Vec3>> color_change = cgogn::get_attribute<Vec3 , Vertex>(m, "color");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m,color_change,v) = BLUE;
			return true;
		});
	}

	void highlight_difference(MESH& m , Attribute<Vec3>* au_position){
		std::shared_ptr<Attribute<Vec3>> color_change = cgogn::get_attribute<Vec3 , Vertex>(m , "color");
		std::shared_ptr<Attribute<Vec3>> pos_au_repos = cgogn::get_attribute<Vec3 , Vertex>(m , "AU00");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			
			if(value<Vec3>(m,au_position,v) != value<Vec3>(m,pos_au_repos,v)){
				value<Vec3>(m,color_change,v) = GREEN;
			}
			else
				value<Vec3>(m,color_change,v) = BLUE;
			return true;
		});
		mesh_provider_->emit_attribute_changed(m,color_change.get());
	}

	void blending(MESH& m , Attribute<Vec3>* attribute_to_blend , float weight){
		std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(m, "position");
		std::shared_ptr<Attribute<Vec3>> color = cgogn::get_attribute<Vec3, Vertex>(m, "color");
		std::shared_ptr<Attribute<Vec3>> repos_position = cgogn::get_attribute<Vec3, Vertex>(m, "AU00");
		Attribute<Vec3>* new_vertex_pos_value = vertex_position.get();
		Vec3 diff_distance_repos = Vec3(0,0,0);

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (value<Vec3>(m,color,v) == GREEN)
			{
				value<Vec3>(m,vertex_position,v) = value<Vec3>(m,repos_position,v);
				diff_distance_repos = value<Vec3>(m,attribute_to_blend,v) - value<Vec3>(m,repos_position,v);
				diff_distance_repos = diff_distance_repos * weight;
				value<Vec3>(m,vertex_position,v) += diff_distance_repos;
			}
			return true;
		});
		mesh_provider_->emit_attribute_changed(m,new_vertex_pos_value);
	}

	void set_distance(MESH &m, Attribute<Vec3>* blendshape_start , Attribute<Vec3>* blendshape_target){
		std::shared_ptr<Attribute<Vec3>> distance = cgogn::get_or_add_attribute<Vec3, Vertex>(m, "distance");
		Attribute<Vec3>* distance_value = distance.get();
		Vec3 diff_distance_repos = Vec3(0,0,0);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			diff_distance_repos = value<Vec3>(m,blendshape_target,v) - value<Vec3>(m,blendshape_start,v);
			value<Vec3>(m,distance_value,v) = diff_distance_repos;
			return true;
		});

	}

	void interpolation(MESH &m, float pas){
		std::shared_ptr<Attribute<Vec3>> distance = cgogn::get_or_add_attribute<Vec3, Vertex>(m, "distance");
		std::shared_ptr<Attribute<Vec3>> position_interpolation = cgogn::get_or_add_attribute<Vec3, Vertex>(m, "position_interpolation");
		Attribute<Vec3>* distance_value = distance.get();
		Attribute<Vec3>* interpolation_value = position_interpolation.get();
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m,interpolation_value,v) += value<Vec3>(m,distance_value,v) * pas ;
			return true;
		});
		mesh_provider_->emit_attribute_changed(m,interpolation_value);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		set_all(directory_,".obj",path_aus_);
		set_all(directory_,".csv",path_csv_);
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
				change_to_selected_au(*selected_mesh_ , selected_vertex_position_.get());

				static const char* current_item_mesh = NULL;
				if (ImGui::BeginCombo("Meshes to blend", current_item_mesh)) // The second parameter is the label previewed before opening the combo.
				{
					
					for (int n = 0; n < pos_aus_.size() ; n++)
					{
						bool is_selected = (current_item_mesh == pos_aus_[n]->name().c_str()); // You can store your selection however you want, outside or inside your objects
						if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected)){
							current_item_mesh = pos_aus_[n]->name().c_str();
							if (!(std::find(std::begin(attribute_to_blend_),std::end(attribute_to_blend_),pos_aus_[n]) != std::end(attribute_to_blend_)))
							{
								attribute_to_blend_.push_back(pos_aus_[n]);
								weights.push_back(1.);
							}
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}

				if (attribute_to_blend_.size() != 0)
				{
					for (int i = 0; i < attribute_to_blend_.size(); i++)
					{
						ImGui::SliderFloat(attribute_to_blend_[i]->name().c_str(), &weight, 0.0, 5.0);
					}
				}
				static bool start = false;
				if (ImGui::Button("Blend"))
				{
					//std::ostringstream new_attribute_name;
					for (int i = 0; i < attribute_to_blend_.size(); i++)
					{
						//new_attribute_name << attribute_to_blend[i]->name().c_str() << 'w' << weights[i] << '+';
						//highlight_difference(*selected_mesh_,attribute_to_blend_[i].get());
						//blending(*selected_mesh_,attribute_to_blend_[i].get(),weights[i]);
					}
					//std::shared_ptr<Attribute<Vec3>> new_attribute = add_attribute<Vec3 , Vertex>(m , new_attribute_name.str().substr(0,new_attribute_name.str().size()-1).c_str());
					
					highlight_difference(*selected_mesh_,attribute_to_blend_[0].get());
					start = true;
					// attribute_to_blend_.clear();
					// weights.clear();
				}

				if (ImGui::Button("Clear"))
				{
					std::shared_ptr<Attribute<Vec3>> repos_position = cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, "AU00");
					blending(*selected_mesh_,repos_position.get(),weight);
					start = false;
					weight = 0;
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
						blending(*selected_mesh_,attribute_to_blend_[0].get(),weight);
						weight += 0.003;
					}
					else
						start = false;
				}
				
				
				

				// ImGui::Separator();

				// static const char* current_item_aus = NULL;
				// if (ImGui::BeginCombo("AUs Difference", current_item_aus)) // The second parameter is the label previewed before opening the combo.
				// {
					
				// 	for (int n = 0; n < pos_aus_.size() ; n++)
				// 	{
				// 		bool is_selected = (current_item_aus == pos_aus_[n]->name().c_str()); // You can store your selection however you want, outside or inside your objects
				// 		if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected)){
				// 			current_item_aus = pos_aus_[n]->name().c_str();
				// 			highlight_difference(*selected_mesh_,pos_aus_[n].get(),RED);
				// 		}
				// 		if (is_selected)
				// 			ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
				// 	}
				// 	ImGui::EndCombo();
				// }

				static const char* current_item_csv = NULL;
				if (ImGui::BeginCombo("Load CSV", current_item_csv)) // The second parameter is the label previewed before opening the combo.
				{
					for (int n = 0; n < path_csv_.size() ; n++)
					{
						bool is_selected = (current_item_csv == path_csv_[n].c_str()); // You can store your selection however you want, outside or inside your objects
						if (ImGui::Selectable(path_csv_[n].c_str(), is_selected)){
							current_item_csv = path_csv_[n].c_str();
							csv_.clear();
							timestamp_csv_.clear();
							csv_parser(path_csv_[n],';');
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
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
						count2 = iter->second.size();
						for (auto &it : csv_)
						{
							if (it.first == "timestamp")
							{
								timestamp_csv_ = it.second;
							}
						}
						for (int i = 0; i < timestamp_csv_.size(); i++)
						{
							std::cout << timestamp_csv_[i] << std::endl;
						}
											
						time_start = ui::App::frame_time_;
						incr = 0.;
						poids_frame = 1.;
						
					}
				}
				
				
				
				if (incr < count2)
				{
					std::vector<int> aus_confirm;
					std::vector<int> aus_confirm_next_frame;
					std::vector<float> weights_frame;
					std::vector<float> weights_next_frame;
					
					for (auto &it : csv_){
						if (ends_with(it.first , "_r"))
						{	
							attributes_csv_.push_back(cgogn::get_attribute<Vec3, Vertex>(*selected_mesh_, it.first.substr(0,it.first.size()-2)));
							if (incr+1 < count2)
							{
								weights_frame.push_back(it.second[incr]);
								weights_next_frame.push_back(it.second[incr+1]);
							}
							else{
								weights_frame.push_back(it.second[incr-1]);
								weights_next_frame.push_back(it.second[incr]);
							}
							if (it.second[incr] > 0.01)
							{
								std::cout << it.first << " poids : " << it.second[incr] << std::endl;
							}
							
							
						}
						if (ends_with(it.first , "_c"))
						{
							if (incr+1 < count2)
							{
								aus_confirm.push_back(it.second[incr]);
								aus_confirm_next_frame.push_back(it.second[incr+1]);
							}
							else{
								aus_confirm.push_back(it.second[incr-1]);
								aus_confirm_next_frame.push_back(it.second[incr]);
							}
						}
					}
					for (int i = 0; i < attributes_csv_.size(); i++)
					{
						if (aus_confirm[i] == 0)
						{
							weights_frame[i] = 0.;
						}
						if (aus_confirm_next_frame[i] == 0)
						{
							weights_next_frame[i] = 0.;
						}
						weights.push_back((weights_frame[i]*(1-poids_frame)) + (weights_next_frame[i]*poids_frame));
						highlight_difference(*selected_mesh_,attributes_csv_[i].get());
						blending(*selected_mesh_,attributes_csv_[i].get(),weights[i]);
					}	
					timer = ui::App::frame_time_ - time_start;

					while((timer > timestamp_csv_[incr]) && (incr < count2))
					{
						incr++;
					}
					if (incr-1 > 0)
					{
						poids_frame = (timer - timestamp_csv_[incr-1])/(timestamp_csv_[incr] - timestamp_csv_[incr-1]);
					}
					std::cout << "timer : " << timer << std::endl;
					std::cout << "poids_frame : " << poids_frame << std::endl;
					std::cout << "timestamp_csv : " << timestamp_csv_[incr] << std::endl;
					
					attributes_csv_.clear();
					weights.clear();
				}

			
				ImGui::Separator();

				ImGui::InputText( "FPS", std::to_string(ui::App::fps()).data() , std::to_string(ui::App::fps()).size());

				ImGui::InputText( "Time since last Frame", std::to_string(ui::App::frame_time_).data() , std::to_string(ui::App::frame_time_).size());

				static const char* current_item_start = NULL;
				if (ImGui::BeginCombo("Interpolation shape start", current_item_start)) // The second parameter is the label previewed before opening the combo.
				{
					
					for (int n = 0; n < pos_aus_.size() ; n++)
					{
						bool is_selected = (current_item_start == pos_aus_[n]->name().c_str()); // You can store your selection however you want, outside or inside your objects
						if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected)){
							current_item_start = pos_aus_[n]->name().c_str();
							au_start = pos_aus_[n].get();
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}

				static const char* current_item_target = NULL;
				if (ImGui::BeginCombo("Interpolation shape target", current_item_target)) // The second parameter is the label previewed before opening the combo.
				{
					
					for (int n = 0; n < pos_aus_.size() ; n++)
					{
						bool is_selected = (current_item_target == pos_aus_[n]->name().c_str()); // You can store your selection however you want, outside or inside your objects
						if (ImGui::Selectable(pos_aus_[n]->name().c_str(), is_selected)){
							current_item_target = pos_aus_[n]->name().c_str();
							au_target = pos_aus_[n].get();
						}
						if (is_selected)
							ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
					}
					ImGui::EndCombo();
				}
				
				ImGui::SliderInt("Nb_frames", &nb_frames, 1, 600);

				if (ImGui::Button("Start Interpolation"))
				{
					if ((current_item_start != NULL) && (current_item_target != NULL))
					{
						set_attribute(*selected_mesh_,au_start,"position_interpolation");
						set_distance(*selected_mesh_,au_start,au_target);
						count = nb_frames;
					}
				}

				if (count != 0)
				{
					interpolation(*selected_mesh_,1./nb_frames);
					count--;
				}
				
			}
		}
	}

private:
	MESH* selected_mesh_;
	View* selected_view_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> pos_aus_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> attribute_to_blend_;
	std::vector<std::shared_ptr<Attribute<Vec3>>> attributes_csv_;
	// std::shared_ptr<Attribute<Vec3>> selected_vertex_normal_;
	MeshProvider<MESH>* mesh_provider_;
	std::vector<std::string> path_aus_;
	std::vector<std::string> path_csv_;
	std::vector<float> weights;
	std::string directory_;
	std::map<std::string,std::vector<float>> csv_;
	std::vector<float> timestamp_csv_;
	Attribute<Vec3>* au_target;
	Attribute<Vec3>* au_start;
	int nb_frames = 1;
	float64 timer = 0.;
	float64 time_start = 0;
	int count = 0;
	int count2 = 0;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_ACTION_UNIT_CHANGE_H_
