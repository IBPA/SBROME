/*
 * main.cpp
 *
 *  Created on: Apr 8, 2012
 *      Author: linh, UC Davis
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <map>

#include "DesignProblem.h"
#include "PartDatabase.h"
#include "Framework.h"

#include "xmlParser.h"

using namespace std;

ValueMatrix transpose (ValueMatrix matrix) {
	ValueMatrix tmp;
	if (!matrix.empty()) {
		int n = matrix.size();
		int m = matrix[0].size();
		tmp.resize(m);
		for (int i = 0; i < m; i++) {
			tmp[i].resize(n);
			 for (int j = 0; j < n; j++)
				 tmp[i][j] = matrix[j][i];
		}
	}
	return tmp;
}

NodeVisualInfo* GetVisualInfo(XMLNode xml_node) {
	return new NodeVisualInfo(atoi(xml_node.getAttribute("x")), atoi(xml_node.getAttribute("y")), atoi(xml_node.getAttribute("width")), atoi(xml_node.getAttribute("height")));
}

void Write2XML(GeneCircuitGraph* gcg, SimulationData desired_behavior, SimulationData simulated_behavior, NodeVisualInfoMap* nvim, string filename) {
	ofstream myfile (filename.c_str());
	if (myfile.is_open()) {
		// WRITE HEADER
		string TEXT1 = "<mxGraphModel grid=\"1\" guides=\"1\" tooltips=\"1\" connect=\"1\" fold=\"1\" page=\"0\" pageScale=\"1\" pageWidth=\"826\" pageHeight=\"1169\">";
		myfile << TEXT1 << "\n";
		TEXT1 = "<root>";
		myfile << TEXT1 << "\n";
		TEXT1 = "<mxCell id=\"0\"/>";
		myfile << TEXT1 << "\n";
		TEXT1 = "<mxCell id=\"1\" parent=\"0\"/>";
		myfile << TEXT1 << "\n";
		int number_of_cells = 0;
		// WRITE NODES
		double graph_border_y = 0;
		for (int node_id = 0; node_id < gcg->NodeCount(); node_id++) {
			number_of_cells++;
			string node_name = gcg->GetNodeAttr(node_id)->component->name;
			XMLNode dummy_node;
			XMLNode xml_node = dummy_node.createXMLTopNode("mxCell");
			xml_node.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_node.addAttribute("value", node_name.c_str());
			string node_style;
			if (gcg->GetNodeAttr(node_id)->component->getCategory() == MOLECULE)
				switch (gcg->GetNodeAttr(node_id)->component->getType()) {
					case m_RNA:
						node_style = "shape=mrna;whiteSpace=wrap;gradientColor=#999999;strokeColor=none";
						break;
					case PROTEIN:
						node_style = "shape=protein;whiteSpace=wrap;gradientColor=#999999;strokeColor=none";
						break;
					case LIGAND:
						node_style = "shape=ligand;whiteSpace=wrap;gradientColor=#999999;strokeColor=none";
						break;
					case PROTEIN_COMPLEX:
						node_style = "shape=protein-protein;whiteSpace=wrap;gradientColor=#999999;strokeColor=none";
						break;
					case LIGAND_PROTEIN_COMPLEX:
						node_style = "shape=ligand-protein;whiteSpace=wrap;gradientColor=#999999;strokeColor=none";
						break;
					case RNA_COMPLEX:
						node_style = "shape=rna-rna;whiteSpace=wrap;gradientColor=#999999;strokeColor=none";
						break;
					case POOL:
						node_style = "round=1";
						break;

				}
			else if (gcg->GetNodeAttr(node_id)->component->getCategory() == LOGIC_GATE) {
				switch (gcg->GetNodeAttr(node_id)->component->getType()) {
					case YES_GATE:
						node_style = "shape=yes_gate";
						break;
					case NOT_GATE:
						node_style = "shape=not_gate";
						break;
					case AND_GATE:
						node_style = "shape=and_gate";
						break;
					case OR_GATE:
						node_style = "shape=or_gate_aux";
						break;
				}
			}
			else if (gcg->GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL) {
				switch (gcg->GetNodeAttr(node_id)->component->getType()) {
					case INPUT:
						node_style = "shape=input";
						break;
					case OUTPUT:
						node_style = "shape=output";
						break;
				}
			}
			if (node_name.length() > 0) {
				double font_size = (*nvim)[node_id]->width/node_name.length();
				node_style += ";fontSize=" + Num2Str(font_size);
			}
			node_style += ";strokeWidth=" + Num2Str((*nvim)[node_id]->line_width);
			xml_node.addAttribute("style", node_style.c_str());
			xml_node.addAttribute("vertex", "1");
			xml_node.addAttribute("parent", "1");
			// Add a geometry child
			XMLNode geometry_node = xml_node.addChild("mxGeometry");
			geometry_node.addAttribute("x", Num2Str((*nvim)[node_id]->pos_x).c_str());
			geometry_node.addAttribute("y", Num2Str((*nvim)[node_id]->pos_y).c_str());
			geometry_node.addAttribute("width", Num2Str((*nvim)[node_id]->width).c_str());
			geometry_node.addAttribute("height", Num2Str((*nvim)[node_id]->height).c_str());
			geometry_node.addAttribute("as","geometry");
			if (graph_border_y < ((*nvim)[node_id]->pos_y + (*nvim)[node_id]->height))
				graph_border_y = (*nvim)[node_id]->pos_y + (*nvim)[node_id]->height;
			char* XMLstr = xml_node.createXMLString();
			myfile << XMLstr << "\n";
			delete []XMLstr;
		}
		// WRITE EDGES
		XMLNode dummy_node;
		for (int src_node_id = 0; src_node_id < gcg->NodeCount(); src_node_id++) {
			for (int dest_node_id = 0; dest_node_id < gcg->NodeCount(); dest_node_id++) {
				if (src_node_id != dest_node_id && gcg->GetEdgeAttr(src_node_id, dest_node_id) != NULL) {
					number_of_cells++;
					XMLNode xml_node = dummy_node.createXMLTopNode("mxCell");
					xml_node.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
					xml_node.addAttribute("value","");
					double line_width = ((*nvim)[src_node_id]->line_width < (*nvim)[dest_node_id]->line_width)?((*nvim)[src_node_id]->line_width):((*nvim)[dest_node_id]->line_width);
					string edge_style = "edgeStyle=elbowEdgeStyle;elbow=horizontal";
					switch (gcg->GetEdgeAttr(src_node_id, dest_node_id)->type) {
						case ACTIVATORY:
							edge_style += ";endArrow=classic";
							break;
						case INHIBITORY:
							edge_style += ";endArrow=oval";
							break;
						case UNKNOWN:
							edge_style += ";endArrow=none";
							break;
					}
					edge_style += ";exitX=1;exitY=0.5;entryX=0;entryY=0.5";
					edge_style += ";strokeWidth=" + Num2Str(line_width) + ";endSize=" + Num2Str(line_width*5);
					xml_node.addAttribute("style",edge_style.c_str());
					xml_node.addAttribute("edge","1");
					xml_node.addAttribute("parent","1");
					xml_node.addAttribute("source",Num2Str(src_node_id + 2).c_str());
					xml_node.addAttribute("target",Num2Str(dest_node_id + 2).c_str());
					// Add a geometry child
					XMLNode geometry_node = xml_node.addChild("mxGeometry");
					geometry_node.addAttribute("width","100");
					geometry_node.addAttribute("height","100");
					geometry_node.addAttribute("as","geometry");
					XMLNode mxpoint_node_1 = geometry_node.addChild("mxPoint");
					mxpoint_node_1.addAttribute("y","100");
					mxpoint_node_1.addAttribute("as","sourcePoint");
					XMLNode mxpoint_node_2 = geometry_node.addChild("mxPoint");
					mxpoint_node_2.addAttribute("x","100");
					mxpoint_node_2.addAttribute("as","targetPoint");
					char* XMLstr = xml_node.createXMLString();
					myfile << XMLstr << "\n";
					delete []XMLstr;
				}
			}
		}
		// WRITE THE SIMULATION DATA
		if (!simulated_behavior.input_values.empty() && !simulated_behavior.output_values.empty()) {
			double data_x = 100;
			double data_y = graph_border_y + 100;
			double data_w = 400;
			double data_h = 200;
			double data_point_h = 10;
			// X_AXIS
			number_of_cells++;
			XMLNode xml_node1 = dummy_node.createXMLTopNode("mxCell");
			xml_node1.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_node1.addAttribute("value", "");
			xml_node1.addAttribute("style", "line;rotation=0;strokeWidth=3");
			xml_node1.addAttribute("vertex", "1");
			xml_node1.addAttribute("parent", "1");
			// Add a geometry child
			XMLNode geometry_node1 = xml_node1.addChild("mxGeometry");
			geometry_node1.addAttribute("x", Num2Str(data_x).c_str());
			geometry_node1.addAttribute("y", Num2Str(data_y + data_h + data_point_h).c_str());
			geometry_node1.addAttribute("width", Num2Str(data_w).c_str());
			geometry_node1.addAttribute("height", "1");
			geometry_node1.addAttribute("as","geometry");
			char* XMLstr = xml_node1.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
			delete []XMLstr;
			// Y_AXIS
			number_of_cells++;
			XMLNode xml_node2 = dummy_node.createXMLTopNode("mxCell");
			xml_node2.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_node2.addAttribute("value", "");
			xml_node2.addAttribute("style", "line;rotation=90;strokeWidth=3");
			xml_node2.addAttribute("vertex", "1");
			xml_node2.addAttribute("parent", "1");
			// Add a geometry child
			XMLNode geometry_node2 = xml_node2.addChild("mxGeometry");
			geometry_node2.addAttribute("x", Num2Str(data_x - data_h/2).c_str());
			geometry_node2.addAttribute("y", Num2Str(data_y + data_point_h).c_str());
			geometry_node2.addAttribute("width", Num2Str(data_h).c_str());
			geometry_node2.addAttribute("height", Num2Str(data_h).c_str());
			geometry_node2.addAttribute("as","geometry");
			XMLstr = xml_node2.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
			delete []XMLstr;
			// DATA POINTS
			double max_concentration = 0;
			for (int row = 0; row < desired_behavior.output_values.size(); row++) {
				if (max_concentration < desired_behavior.output_values[row][0])
					max_concentration = desired_behavior.output_values[row][0];
				if (max_concentration < simulated_behavior.output_values[row][0])
					max_concentration = simulated_behavior.output_values[row][0];
			}
			// X_LABEL
			for (int col = 0; col < desired_behavior.input_values[0].size(); col++) {
				number_of_cells++;
				XMLNode xml_x_label = dummy_node.createXMLTopNode("mxCell");
				xml_x_label.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
				string tmp = "[" + gcg->GetNodeAttr(desired_behavior.inputs[col])->component->name + "]";
				xml_x_label.addAttribute("value", tmp.c_str());
				xml_x_label.addAttribute("style", "text");
				xml_x_label.addAttribute("vertex", "1");
				xml_x_label.addAttribute("parent", "1");
				// Add a geometry child
				XMLNode geometry_x_label = xml_x_label.addChild("mxGeometry");
				geometry_x_label.addAttribute("x", Num2Str(data_x + 5).c_str());
				geometry_x_label.addAttribute("y", Num2Str(data_y + data_h + data_point_h + col*20 + 20).c_str());
				geometry_x_label.addAttribute("width", "10");
				geometry_x_label.addAttribute("height", "10");
				geometry_x_label.addAttribute("as","geometry");
				XMLstr = xml_x_label.createXMLString();
				cout << XMLstr << endl;
				myfile << XMLstr << "\n";
				delete []XMLstr;
			}
			for (int row = 0; row < desired_behavior.output_values.size(); row++) {
				// X_STICK
				for (int col = 0; col < desired_behavior.input_values[0].size(); col++) {
					number_of_cells++;
					XMLNode xml_x_label = dummy_node.createXMLTopNode("mxCell");
					xml_x_label.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
					xml_x_label.addAttribute("value", Num2Str(desired_behavior.input_values[row][col]).c_str());
					xml_x_label.addAttribute("style", "text");
					xml_x_label.addAttribute("vertex", "1");
					xml_x_label.addAttribute("parent", "1");
					// Add a geometry child
					XMLNode geometry_x_label = xml_x_label.addChild("mxGeometry");
					geometry_x_label.addAttribute("x", Num2Str(data_w*(row + 1)/(desired_behavior.output_values.size() + 1) + data_x).c_str());
					geometry_x_label.addAttribute("y", Num2Str(data_y + data_h + data_point_h + col*20 + 20).c_str());
					geometry_x_label.addAttribute("width", "10");
					geometry_x_label.addAttribute("height", "10");
					geometry_x_label.addAttribute("as","geometry");
					XMLstr = xml_x_label.createXMLString();
					cout << XMLstr << endl;
					myfile << XMLstr << "\n";
					delete []XMLstr;
				}
				// DESIRED VALUE
				number_of_cells++;
				XMLNode xml_node1 = dummy_node.createXMLTopNode("mxCell");
				xml_node1.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
				xml_node1.addAttribute("value", "");
				xml_node1.addAttribute("style", "ellipse;whiteSpace=wrap");
				xml_node1.addAttribute("vertex", "1");
				xml_node1.addAttribute("parent", "1");
				// Add a geometry child
				XMLNode geometry_node1 = xml_node1.addChild("mxGeometry");
				geometry_node1.addAttribute("x", Num2Str(data_w*(row + 1)/(desired_behavior.output_values.size() + 1) + data_x).c_str());
				geometry_node1.addAttribute("y", Num2Str(0.9*data_h*desired_behavior.output_values[row][0]/max_concentration + data_y).c_str());
				geometry_node1.addAttribute("width", "10");
				geometry_node1.addAttribute("height", "10");
				geometry_node1.addAttribute("as","geometry");
				XMLstr = xml_node1.createXMLString();
				cout << XMLstr << endl;
				myfile << XMLstr << "\n";
				delete []XMLstr;
				// SIMULATED VALUE
				number_of_cells++;
				XMLNode xml_node2 = dummy_node.createXMLTopNode("mxCell");
				xml_node2.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
				xml_node2.addAttribute("value", "");
				xml_node2.addAttribute("style", "rhombus;whiteSpace=wrap");
				xml_node2.addAttribute("vertex", "1");
				xml_node2.addAttribute("parent", "1");
				// Add a geometry child
				XMLNode geometry_node2 = xml_node2.addChild("mxGeometry");
				geometry_node2.addAttribute("x", Num2Str(data_w*(row + 1)/(desired_behavior.output_values.size() + 1) + data_x).c_str());
				geometry_node2.addAttribute("y", Num2Str(0.9*data_h*simulated_behavior.output_values[row][0]/max_concentration + data_y).c_str());
				geometry_node2.addAttribute("width", "10");
				geometry_node2.addAttribute("height", "10");
				geometry_node2.addAttribute("as","geometry");
				XMLstr = xml_node2.createXMLString();
				cout << XMLstr << endl;
				myfile << XMLstr << "\n";
				delete []XMLstr;
			}
			// Y_LABEL
			number_of_cells++;
			XMLNode xml_y_label = dummy_node.createXMLTopNode("mxCell");
			xml_y_label.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			string tmp = "[" + gcg->GetNodeAttr(desired_behavior.outputs[0])->component->name + "]";
			xml_y_label.addAttribute("value", tmp.c_str());
			xml_y_label.addAttribute("style", "text");
			xml_y_label.addAttribute("vertex", "1");
			xml_y_label.addAttribute("parent", "1");
			// Add a geometry child
			XMLNode geometry_y_label = xml_y_label.addChild("mxGeometry");
			geometry_y_label.addAttribute("x", Num2Str(data_x - 40).c_str());
			geometry_y_label.addAttribute("y", Num2Str(data_y + data_point_h).c_str());
			geometry_y_label.addAttribute("width", "10");
			geometry_y_label.addAttribute("height", "10");
			geometry_y_label.addAttribute("as","geometry");
			XMLstr = xml_y_label.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
			delete []XMLstr;
			// LEGEND
			// DESIRED LEGEND
			number_of_cells++;
			XMLNode xml_desired_legend = dummy_node.createXMLTopNode("mxCell");
			xml_desired_legend.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_desired_legend.addAttribute("value", "Desired");
			xml_desired_legend.addAttribute("style", "text");
			xml_desired_legend.addAttribute("vertex", "1");
			xml_desired_legend.addAttribute("parent", "1");
			// Add a geometry child
			XMLNode geometry_desired_legend = xml_desired_legend.addChild("mxGeometry");
			geometry_desired_legend.addAttribute("x", Num2Str(data_x + data_w + 20).c_str());
			geometry_desired_legend.addAttribute("y", Num2Str(data_y - 20).c_str());
			geometry_desired_legend.addAttribute("width", "10");
			geometry_desired_legend.addAttribute("height", "10");
			geometry_desired_legend.addAttribute("as","geometry");
			XMLstr = xml_desired_legend.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
			delete []XMLstr;
			// DESIRED ICON
			number_of_cells++;
			xml_desired_legend = dummy_node.createXMLTopNode("mxCell");
			xml_desired_legend.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_desired_legend.addAttribute("value", "");
			xml_desired_legend.addAttribute("style", "ellipse;whiteSpace=wrap");
			xml_desired_legend.addAttribute("vertex", "1");
			xml_desired_legend.addAttribute("parent", "1");
			// Add a geometry child
			geometry_desired_legend = xml_desired_legend.addChild("mxGeometry");
			geometry_desired_legend.addAttribute("x", Num2Str(data_x + data_w).c_str());
			geometry_desired_legend.addAttribute("y", Num2Str(data_y + data_point_h - 20).c_str());
			geometry_desired_legend.addAttribute("width", "10");
			geometry_desired_legend.addAttribute("height", "10");
			geometry_desired_legend.addAttribute("as","geometry");
			XMLstr = xml_desired_legend.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
			delete []XMLstr;
			// SIMULATION LEGEND
			number_of_cells++;
			xml_desired_legend = dummy_node.createXMLTopNode("mxCell");
			xml_desired_legend.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_desired_legend.addAttribute("value", "Simulation");
			xml_desired_legend.addAttribute("style", "text");
			xml_desired_legend.addAttribute("vertex", "1");
			xml_desired_legend.addAttribute("parent", "1");
			// Add a geometry child
			geometry_desired_legend = xml_desired_legend.addChild("mxGeometry");
			geometry_desired_legend.addAttribute("x", Num2Str(data_x + data_w + 20).c_str());
			geometry_desired_legend.addAttribute("y", Num2Str(data_y).c_str());
			geometry_desired_legend.addAttribute("width", "10");
			geometry_desired_legend.addAttribute("height", "10");
			geometry_desired_legend.addAttribute("as","geometry");
			XMLstr = xml_desired_legend.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
			delete []XMLstr;
			// SIMULATION ICON
			number_of_cells++;
			xml_desired_legend = dummy_node.createXMLTopNode("mxCell");
			xml_desired_legend.addAttribute("id", Num2Str(number_of_cells + 1).c_str());
			xml_desired_legend.addAttribute("value", "");
			xml_desired_legend.addAttribute("style", "rhombus;whiteSpace=wrap");
			xml_desired_legend.addAttribute("vertex", "1");
			xml_desired_legend.addAttribute("parent", "1");
			// Add a geometry child
			geometry_desired_legend = xml_desired_legend.addChild("mxGeometry");
			geometry_desired_legend.addAttribute("x", Num2Str(data_x + data_w).c_str());
			geometry_desired_legend.addAttribute("y", Num2Str(data_y + data_point_h).c_str());
			geometry_desired_legend.addAttribute("width", "10");
			geometry_desired_legend.addAttribute("height", "10");
			geometry_desired_legend.addAttribute("as","geometry");
			XMLstr = xml_desired_legend.createXMLString();
			cout << XMLstr << endl;
			myfile << XMLstr << "\n";
		}
		// WRITE FOOTER
		myfile << "</root>" << "\n" << "</mxGraphModel>";
    }
	myfile.close();
}

// void Run(string database_file, string input_file)
void Run(string filename) {
	PartDatabase* pdb = new PartDatabase("Database/CAD");
	vector<Module*> *module_lib = pdb->loadModuleLibrary();
	DesignProblem dp;
	ARGEdit editor;
	// PARSE THE INPUT FILE
	XMLNode xMainNode=XMLNode::openFileHelper(filename.c_str());
	int number_of_cells = xMainNode.getChildNode(0).nElement();
	int number_of_nodes = 0;
	int number_of_inputs = 0;
	map<int,int> cell_to_node;
	map<int,NodeVisualInfo*> node_id_to_visual_info;
	ValueMatrix input_values_t, output_values_t;	// transpose
	IdList inputs, outputs;
	ValueList vl_tmp;
	// BUILD THE INPUT GRAPH
	// Add nodes
	for (int i = 2; i < number_of_cells; i++) {
		int cell_id = atoi(xMainNode.getChildNode(0).getChildNode(i).getAttribute("id"));
		XMLNode current_node;
		if(strcmp(xMainNode.getChildNode(0).getChildNode(i).getName(), "mxCell") == STR_EQ)
			current_node = xMainNode.getChildNode(0).getChildNode(i);
		else if (strcmp(xMainNode.getChildNode(0).getChildNode(i).getName(), "UserObject") == STR_EQ) {	// input or output
			vl_tmp.clear();
			string values = xMainNode.getChildNode(0).getChildNode(i).getAttribute("link");
			int pre_pos = -1;
			for (int tmp = 0; tmp < values.length(); tmp++) {
				if (values[tmp] != ' ' && pre_pos < 0)
					pre_pos = tmp;
				if (values[tmp] == ' ' && pre_pos >= 0) {
					vl_tmp.push_back(atof(values.substr(pre_pos, tmp - pre_pos).c_str()));
					pre_pos = -1;
				}
				if (tmp == values.length() - 1 && pre_pos >= 0)
					vl_tmp.push_back(atof(values.substr(pre_pos, tmp - pre_pos + 1).c_str()));
			}
			current_node = xMainNode.getChildNode(0).getChildNode(i).getChildNode(0);
		}
		string style_str = current_node.getAttribute("style");
		//cout << style_str << endl;
		string l_str = "shape=ligand;", p_str = "shape=protein;", m_str = "shape=mrna;";
		string l_p_str = "shape=ligand-protein;", p_p_str = "shape=protein-protein;", m_m_str = "shape=mrna-mrna;";
		if (style_str.compare("shape=not_gate") == STR_EQ) {						// 	NOT gate
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("NOT_GATE"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.compare("shape=yes_gate") == STR_EQ) {					// 	NOT gate
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("YES_GATE"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.compare("shape=or_gate") == STR_EQ) {					//	OR gate
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("OR_GATE"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.compare("shape=and_gate") == STR_EQ) {					//	AND gate
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("AND_GATE"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.compare("shape=input") == STR_EQ) {						//	input
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("INPUT"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
			inputs.push_back(cell_to_node[cell_id]);
			input_values_t.push_back(vl_tmp);
		}
		else if (style_str.compare("shape=output") == STR_EQ) {						//	output
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("OUTPUT"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
			outputs.push_back(cell_to_node[cell_id]);
			output_values_t.push_back(vl_tmp);
		}
		else if (style_str.substr(0,l_str.length()).compare(l_str) == STR_EQ) {		// Ligand
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("Ligand"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.substr(0,p_str.length()).compare(p_str) == STR_EQ) {		// Protein
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("Protein"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.substr(0,m_str.length()).compare(m_str) == STR_EQ) {		// mRNA
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("mRNA"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.substr(0,l_p_str.length()).compare(l_p_str) == STR_EQ) {	// Ligand-Protein
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("LigandProteinComplex"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.substr(0,p_p_str.length()).compare(p_p_str) == STR_EQ) {	// Protein-Protein
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("ProteinComplex"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
		else if (style_str.substr(0,m_m_str.length()).compare(m_m_str) == STR_EQ) {	// RNA-RNA
			cell_to_node[cell_id] = editor.InsertNode(CreateBioNetNode("RNAComplex"));
			node_id_to_visual_info[cell_to_node[cell_id]] = GetVisualInfo(current_node.getChildNode(0));
		}
	}
	// Add edges
	for (int i = 2; i < number_of_cells; i++) {
		XMLNode current_node;
		if(strcmp(xMainNode.getChildNode(0).getChildNode(i).getName(), "mxCell") == STR_EQ)
			current_node = xMainNode.getChildNode(0).getChildNode(i);
		else if (strcmp(xMainNode.getChildNode(0).getChildNode(i).getName(), "UserObject") == STR_EQ)
			current_node = xMainNode.getChildNode(0).getChildNode(i).getChildNode(0);
		if (strncmp(current_node.getAttribute("style") , "edgeStyle", strlen("edgeStyle")) == STR_EQ) {
			int src_id = atoi(xMainNode.getChildNode(0).getChildNode(i).getAttribute("source"));
			int dest_id = atoi(xMainNode.getChildNode(0).getChildNode(i).getAttribute("target"));
			string edge_str = xMainNode.getChildNode(0).getChildNode(i).getAttribute("style");
			int edge_type = UNKNOWN;
			if (edge_str.find("endArrow=classic") != string::npos)
				edge_type = ACTIVATORY;
			else if (edge_str.find("endArrow=oval") != string::npos)
				edge_type= INHIBITORY;
			editor.InsertEdge(cell_to_node[src_id], cell_to_node[dest_id], new BioNetEdge(edge_type));
		}
	}
	// EXECUTE
	dp.gcg = *(new GeneCircuitGraph(&editor));
	dp.desired_behavior.inputs = inputs;
	dp.desired_behavior.outputs = outputs;
	dp.desired_behavior.input_values = transpose(input_values_t);
	dp.desired_behavior.output_values = transpose(output_values_t);
	//dp.gcg.TestPrint(&std::cout);
	dp.TestPrint(&std::cout);
	Framework framework(pdb,dp);
	vector<std::pair<GeneCircuitGraph*, SimulationData> >* solution_list = framework.Optimize(SI);
	if (!solution_list->empty()) {
		//solution_list->at(0).second.TestPrint(&std::cout);
		GeneCircuitGraph* result_gcg = solution_list->at(0).first;
		// WRITE THE OUTPUT FILE
		// Expand the graph
		GeneCircuitGraph* final_gcg = Layout(result_gcg, pdb, &node_id_to_visual_info);
		//final_gcg->TestPrint(&std::cout);
		// Write into an xml file
		Write2XML(final_gcg, dp.desired_behavior, solution_list->at(0).second, &node_id_to_visual_info, "Database/XML_ORDERED_SCALED_FOR_OUTPUT.xml");
	}
	else
		cout << "No solution!" << endl;
}

double Test(string filename, int module_size) {
	//cout << "##############" << filename << "##############" << endl;
	clock_t start_time = clock();
	PartDatabase* pdb;
	switch (module_size) {
		case 100:
			pdb = new PartDatabase("Database/CAD_100.sqlite");
			break;
		case 50:
			pdb = new PartDatabase("Database/CAD_50.sqlite");
			break;
		case 25:
			pdb = new PartDatabase("Database/CAD_25.sqlite");
			break;
	}
	DesignProblem dp(filename);
	Framework framework(pdb,dp);
	framework.Optimize(GA);
	clock_t end_time = clock();
	return ((double)(end_time - start_time))/((double)CLOCKS_PER_SEC);
}

void execute(string input_filename) {

}

int main() {
	//Run("Database/Input.xml");
	//Run("Benchmarks/Simple/NOT_AND_full.xml");
	//Run("Benchmarks/Simple/NOT_full_2.xml");
	//Run("Benchmarks/Simple/AND.xml");
	//Run("Benchmarks/Simple/AND_full_1.xml");
	//Run("Benchmarks/2-to-1-multiplexer-1.xml");
	Run("Benchmarks/Simple/NOT_full_11.xml");
	//Run("Benchmarks/Simple/NOT_full_12.xml");
	cout << "FINISH!" << endl;
}
