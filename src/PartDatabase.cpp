/*
 * PartDatabase.cpp
 *
 *  Created on: Aug 8, 2012
 *      Author: linh
 */

#include <iostream>
#include <map>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PartDatabase.h"

static int callback_count(void *NotUsed, int argc, char **argv, char **azColName){
	(*(int*)NotUsed) = atoi(argv[0]);
	return 0;
}
static int callback_string(void *NotUsed, int argc, char **argv, char **azColName){
	vector<string> tmp(argc);
	for (int i = 0; i < argc; i++) {
		string str_tmp(argv[i]);
		tmp[i] = str_tmp;
	}
	((StringMatrix*)NotUsed)->push_back(tmp);
	return 0;
}
PartDatabase::PartDatabase() {
	db = NULL;
}

PartDatabase::PartDatabase(std::string db_filename) {
	int rc = sqlite3_open(db_filename.c_str(), &db);
	if (rc) {
		std::cout << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
	}
}

PartDatabase::~PartDatabase() {
	sqlite3_close(db);
}

int PartDatabase::getNumberOfProteins() const {
	return number_of_proteins;
}

double PartDatabase::getDegradationRate(string protein_name) const {
	QueryResultInfo qri = Init();
	qri.query = "SELECT degradation_rate FROM Protein";
	qri.query += " WHERE name = '" + protein_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &qri.string_matrix, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (qri.string_matrix.empty())
		cout << "There is no protein " << protein_name << endl;
	return atof(qri.string_matrix[0][0].c_str());
}

int PartDatabase::getNumberOfMutants(string promoter_name) const {
	QueryResultInfo qri = Init();
	qri.query = "SELECT COUNT(*) FROM Promoter_Mutant";
	qri.query += " WHERE promoter_name = '" + promoter_name + "'";
	int count_tmp;
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_count, &count_tmp, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	return count_tmp;
}

double PartDatabase::getPromoterBasal(string promoter_name, int mutant_index) const {
	QueryResultInfo qri = Init();
	qri.query = "SELECT basal FROM Promoter_Mutant";
	qri.query += " WHERE promoter_name = '" + promoter_name + "' AND mutant_id = " + Num2Str(mutant_index);
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &qri.string_matrix, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (qri.string_matrix.empty())
		cout << "There is no mutant of " << promoter_name << " with mutant index is " << mutant_index << endl;
	return atof(qri.string_matrix[0][0].c_str());
}
double PartDatabase::getPromoterStrength(string promoter_name, int mutant_index) const {
	QueryResultInfo qri = Init();
	qri.query = "SELECT strength FROM Promoter_Mutant";
	qri.query += " WHERE promoter_name = '" + promoter_name + "' AND mutant_id = " + Num2Str(mutant_index);
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &qri.string_matrix, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (qri.string_matrix.empty())
		cout << "There is no mutant of " << promoter_name << " with mutant index is " << mutant_index << endl;
	return atof(qri.string_matrix[0][0].c_str());
}
TranscriptionRegulation PartDatabase::getPromoterRegulation(string regulator_name, string promoter_name) const {
	TranscriptionRegulation tr;
	QueryResultInfo qri = Init();
	qri.query = "SELECT beta,eta,type FROM Transcription_Regulation";
	qri.query += " WHERE regulator_name = '" + regulator_name + "' AND promoter_name = '" + promoter_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &qri.string_matrix, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (qri.string_matrix.empty())
		cout << "There is no interaction between " << regulator_name << " and " << promoter_name << endl;
	tr.binding_affinity = atof(qri.string_matrix[0][0].c_str());
	tr.cooperativity = atof(qri.string_matrix[0][1].c_str());
	tr.regulator_type = (qri.string_matrix[0][2].compare("ACT") == STR_EQ) ? ACTIVATORY : INHIBITORY;
	return tr;
}
LigandRegulation PartDatabase::getLigandRegulation(string ligand_name, string protein_name) const {
	LigandRegulation lr;
	QueryResultInfo qri = Init();
	qri.query = "SELECT theta,eta FROM Ligand_Protein_Interaction";
	qri.query += " WHERE ligand_name = '" + ligand_name + "' AND protein_name = '" + protein_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &qri.string_matrix, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	lr.dissociation = atof(qri.string_matrix[0][0].c_str());
	lr.cooperativity = atof(qri.string_matrix[0][1].c_str());
	return lr;
}

int IOStr2Num(string io_str) {
	if (io_str.compare("") == STR_EQ)
		return INTERMEDIATE;
	else if (io_str.compare("Input") == STR_EQ)
		return INPUT;
	else
		return OUTPUT;
}

vector<Module*>* PartDatabase::loadModuleLibrary() {
	StringMatrix module_info_list;
	QueryResultInfo qri = Init();
	qri.query = "SELECT module_name, motif_name, ref FROM module_info WHERE type = 'FULL'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_info_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	vector<Module*>* module_list = new vector<Module*>(module_info_list.size());

	for (unsigned int module_id = 0; module_id < module_list->size(); module_id++) {
		ARGEdit* module_editor = new ARGEdit;
		IdList inputs, outputs;
		inputs.clear();
		outputs.clear();
		// Read node name
		StringMatrix node_name_list;
		Reset(qri);
		qri.query = "SELECT node_id, node_name FROM module_node";
		qri.query += " WHERE module_name = '" + module_info_list[module_id][0] + "'";
		qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &node_name_list, &qri.zErrMsg);
		if(qri.rc != SQLITE_OK) {
			cout << qri.query << endl;
			std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
			sqlite3_free(qri.zErrMsg);
		}
		map<string, string> node_id_to_node_name;
		for (unsigned int node_id = 0; node_id < node_name_list.size(); node_id++)
			node_id_to_node_name[node_name_list[node_id][0]] = node_name_list[node_id][1];
		// Read nodes
		StringMatrix module_node_list;
		Reset(qri);
		qri.query = "SELECT node_id, node_type, node_role FROM motif_node";
		qri.query += " WHERE motif_name = '" + module_info_list[module_id][1] + "'";
		qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_node_list, &qri.zErrMsg);
		if(qri.rc != SQLITE_OK) {
			cout << qri.query << endl;
			std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
			sqlite3_free(qri.zErrMsg);
		}
		// Read all information
		std::map<string, int> node_id_to_bionode_id;
		for (unsigned int node_id = 0; node_id < module_node_list.size(); node_id++) {
			BioNetNode* bnn = CreateBioNetNode(module_node_list[node_id][1], node_id_to_node_name[module_node_list[node_id][0]]);
			node_id_to_bionode_id[module_node_list[node_id][0]] = module_editor->InsertNode(bnn);
			if (module_node_list[node_id][2].compare("Input") == STR_EQ)
				inputs.push_back(node_id_to_bionode_id[module_node_list[node_id][0]]);
			else if (module_node_list[node_id][2].compare("Output") == STR_EQ)
				outputs.push_back(node_id_to_bionode_id[module_node_list[node_id][0]]);
			else if (module_node_list[node_id][2].compare("Input-Output") == STR_EQ) {
				inputs.push_back(node_id_to_bionode_id[module_node_list[node_id][0]]);
				outputs.push_back(node_id_to_bionode_id[module_node_list[node_id][0]]);
			}
		}
		// Read edges
		StringMatrix module_edge_list;
		Reset(qri);
		qri.query = "SELECT src_node_id, dest_node_id, edge_type FROM motif_edge";
		qri.query += " WHERE motif_name = '" + module_info_list[module_id][1] + "'";
		qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_edge_list, &qri.zErrMsg);
		if(qri.rc != SQLITE_OK) {
			cout << qri.query << endl;
			std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
			sqlite3_free(qri.zErrMsg);
		}
		for (unsigned int edge_id = 0; edge_id < module_edge_list.size(); edge_id++)
			module_editor->InsertEdge(node_id_to_bionode_id[module_edge_list[edge_id][0]],node_id_to_bionode_id[module_edge_list[edge_id][1]], new BioNetEdge(module_edge_list[edge_id][2]));
		// Process for the exceptional cases, it should be handled by a language rather than a graph
		if (module_info_list[module_id][1].compare("input_ligand") == STR_EQ) { // For the case where the module only contains one node
			int input_node_id = module_editor->InsertNode(CreateBioNetNode("INPUT"));
			((InputSignal*)((BioNetNode*)module_editor->GetNodeAttr(input_node_id))->component)->physical_instance = ((BioNetNode*)module_editor->GetNodeAttr(outputs[0]))->component->clone();
			((BioNetNode*)module_editor->GetNodeAttr(input_node_id))->component->name = ((BioNetNode*)module_editor->GetNodeAttr(outputs[0]))->component->name;
			inputs.push_back(input_node_id);
			module_editor->InsertEdge(input_node_id,outputs[0],new BioNetEdge);
		}
		else if (module_info_list[module_id][1].compare("protein_output") == STR_EQ) {
			int output_node_id = module_editor->InsertNode(CreateBioNetNode("OUTPUT"));
			((OutputSignal*)((BioNetNode*)module_editor->GetNodeAttr(output_node_id))->component)->physical_instance = ((BioNetNode*)module_editor->GetNodeAttr(inputs[0]))->component->clone();
			((BioNetNode*)module_editor->GetNodeAttr(output_node_id))->component->name = ((BioNetNode*)module_editor->GetNodeAttr(inputs[0]))->component->name;
			outputs.push_back(output_node_id);
			module_editor->InsertEdge(inputs[0],output_node_id,new BioNetEdge);
		}
		// Create a new module
		module_list->at(module_id) = new Module(module_editor, inputs, outputs);
		module_list->at(module_id)->name = module_info_list[module_id][0];
		// TEST HERE
		//module_list->at(module_id)->TestPrint(&std::cout);
	}
	vector<Module*>* module_list_tmp = GenerateModuleLibrary("Database/ML_June_2013.txt");
	//for (int i = 0 ; i < module_list_tmp->size(); i++)
	//	module_list->push_back(module_list_tmp->at(i));
	delete module_list_tmp;
	return module_list;
}

vector<Module*>* PartDatabase::GenerateModuleLibrary(string filename) {
	vector<Module*>* module_list = new vector<Module*>;
	StringMatrix* ligand_production_list = getAllInteractionSpecies("Ligand");
	StringMatrix* complex_production_list = getAllInteractionSpecies("LigandProteinComplex");
	AppendStringMatrix(complex_production_list, getAllInteractionSpecies("RNAComplex"));
	AppendStringMatrix(complex_production_list, getAllInteractionSpecies("ProteinProteinComplex"));
	ifstream infile;
	infile.open(filename.c_str());
	while(!infile.eof()) {
		string module_name, part_list, inputs, outputs, author, title, journal, year;
		infile >> module_name;
		if (module_name.empty())
			break;	// ???
		infile >> part_list;
		vector<string>* species_name_list = split_string(part_list,'-');
		//cout << part_list << endl;
		infile >> inputs;
		vector<string>* input_list = split_string(inputs, ',');
		infile >> outputs;
		vector<string>* output_list = split_string(outputs, ',');
		// Analyze the part list here
		ARGEdit* module_editor = new ARGEdit;
		int current_promoter_node_id = UNKNOWN;
		map<string,int> available_species;
		IdList direct_product_node_id_list;
		int* node_id_map = new int[species_name_list->size()];	// map from species id to node id in the graph
		for (int species_id = 0; species_id < species_name_list->size(); species_id++) {
			string species_type = getSpeciesType(species_name_list->at(species_id));
			int node_id_tmp =  module_editor->InsertNode(CreateBioNetNode(species_type, species_name_list->at(species_id)));
			available_species[species_name_list->at(species_id)] = node_id_map[species_id] = node_id_tmp;
			direct_product_node_id_list.push_back(node_id_tmp);
			if (species_type.compare("mRNA") == STR_EQ || species_type.compare("mRNA(t)") == STR_EQ)
				current_promoter_node_id = node_id_map[species_id];
			else
				module_editor->InsertEdge(current_promoter_node_id, node_id_map[species_id], new BioNetEdge("ACT"));
		}
		// Add input nodes
		IdList input_node_id_list;
		if (input_list->at(0).compare("EMPTY") != STR_EQ)
			for (int i = 0; i < input_list->size(); i++)
				if (available_species.find(input_list->at(i)) == available_species.end()) {
					int node_id_tmp = module_editor->InsertNode(CreateBioNetNode(getSpeciesType(input_list->at(i)), input_list->at(i)));
					available_species[input_list->at(i)] = node_id_tmp;
					input_node_id_list.push_back(node_id_tmp);
				}
				else
					input_node_id_list.push_back(available_species.find(input_list->at(i))->second);
		// Add ligand nodes which are produced from proteins
		for (int i = 0; i < ligand_production_list->size(); i++)
			if (available_species.find(ligand_production_list->at(i)[0]) == available_species.end()) {
				int j;
				for (j = 1; j < ligand_production_list->at(i).size(); j++)
					if (available_species.find(ligand_production_list->at(i)[j]) == available_species.end())
						break;
				if (j == ligand_production_list->at(i).size()) {
					int node_id_tmp = module_editor->InsertNode(CreateBioNetNode(getSpeciesType(ligand_production_list->at(i)[0]), ligand_production_list->at(i)[0]));
					available_species[ligand_production_list->at(i)[0]] = node_id_tmp;
				}
			}
		// Add complex nodes which are produced from proteins
		IdList complex_node_id_list;
		for (int i = 0; i < complex_production_list->size(); i++) {
			if (available_species.find(complex_production_list->at(i)[0]) == available_species.end()) {
				//cout << complex_production_list->at(i)[0] << endl;
				int j;
				for (j = 1; j < complex_production_list->at(i).size(); j++)
					if (available_species.find(complex_production_list->at(i)[j]) == available_species.end()) {
						//cout << complex_production_list->at(i)[j] << endl;
						break;
					}
				if (j == complex_production_list->at(i).size()) {
					int node_id_tmp = module_editor->InsertNode(CreateBioNetNode(getSpeciesType(complex_production_list->at(i)[0]), complex_production_list->at(i)[0]));
					available_species[complex_production_list->at(i)[0]] = node_id_tmp;
					complex_node_id_list.push_back(node_id_tmp);
					complex_node_id_list.push_back(node_id_tmp);
				}
			}
		}
		// Add output nodes
		IdList output_node_id_list;
		if (output_list->at(0).compare("EMPTY") != STR_EQ)
			for (int i = 0; i < output_list->size(); i++)
				if (available_species.find(output_list->at(i)) != available_species.end())
					output_node_id_list.push_back(available_species.find(output_list->at(i))->second);
				else
					cout << "SBROME: Error! the output " << output_list->at(i) << " is irrelavant with module " << module_name << endl;
		// Add edges
		for (int src_node_id = 0; src_node_id < module_editor->NodeCount(); src_node_id++)
			for (int dest_node_id = 0; dest_node_id < module_editor->NodeCount(); dest_node_id++) {
				string edge_type = getSpeciesInteractionType(((BioNetNode*)module_editor->GetNodeAttr(src_node_id))->component->name,((BioNetNode*)module_editor->GetNodeAttr(dest_node_id))->component->name);
				if (edge_type.compare("NONE") != STR_EQ)
					module_editor->InsertEdge(src_node_id, dest_node_id, new BioNetEdge(edge_type));
			}
		// End analysis
		author = ReadSentence(&infile);
		title = ReadSentence(&infile);
		journal = ReadSentence(&infile);
		infile >> year;
		string ref = author + ". " + title + ". " + journal + ", " + year;
		Module* original_module = new Module(module_editor, input_node_id_list, output_node_id_list);
		//cout << "+++++++++++++++++" << endl;
		//original_module->TestPrint(&std::cout);
		//cout << "+++++++++++++++++" << endl;
		// Generate all fragment modules here
		// Mark ligand nodes
		bool* is_direct_product_node = new bool[original_module->NodeCount()];
		for (int i = 0; i < original_module->NodeCount(); i++)
			is_direct_product_node[i] = false;
		for (int i = 0; i < direct_product_node_id_list.size(); i++)
			is_direct_product_node[i] = true;
		// Mark complex nodes
		bool* is_complex_node = new bool[original_module->NodeCount()];
		for (int i = 0; i < original_module->NodeCount(); i++)
			is_complex_node[i] = false;
		for (int i = 0; i < complex_node_id_list.size(); i++)
			is_complex_node[i] = true;
		// Generate all fragments
		int minimum_length = 2;	// threshold value for the fragment size
		for (int left_i = 0; left_i < species_name_list->size() - minimum_length + 1; left_i++)
			for (int right_i = left_i + minimum_length - 1; right_i < species_name_list->size(); right_i++) {
				// clone here
				ARGEdit* tmp_editor = (ARGEdit*) original_module->ConvertToLoader();
				bool* is_removed = new bool [tmp_editor->NodeCount()];
				// Determine gene product nodes
				for (int i = 0; i < tmp_editor->NodeCount(); i++)
					is_removed[i] = true;
				for (int i = left_i; i <= right_i; i++)
					is_removed[node_id_map[i]] = false;
				// Determine complex nodes
				bool is_continue = true;
				while (is_continue) {
					is_continue = false;
					for (int i = 0; i < original_module->NodeCount(); i++) {
						if (is_removed[i] && original_module->InEdgeCount(i) && !is_direct_product_node[i]) {
							bool is_produced_completely = true;
							bool is_produced_partially = false;
							bool is_used = false;
							for (int j = 0; j < original_module->NodeCount(); j++)
								if (i != j && original_module->GetEdgeAttr(j,i) != NULL)
									if (is_removed[j]) {
										is_produced_completely = false;
										break;
									}
									else
										is_produced_partially = true;
							for (int j = 0; j < original_module->NodeCount(); j++)
								if (i != j && original_module->GetEdgeAttr(i,j) != NULL && !is_removed[j]) {
									is_used = true;
									break;
								}
							if (is_produced_completely || (is_produced_partially && is_used)) {
								is_removed[i] = false;
								is_continue = true;
							}
						}
					}
				}
				IdList input_node_id_list_tmp, output_node_id_list_tmp;
				// Determine the inputs
				for (int i = 0 ; i < original_module->NodeCount(); i++) {
					if (is_removed[i]) {
						int out_edge_num = original_module->OutEdgeCount(i);
						for	(int j = 0; j < out_edge_num; j++) {
							if (!is_removed[original_module->GetOutEdge(i,j)]) {
								input_node_id_list_tmp.push_back(i);
								break;
							}
						}
					}
				}
				for (int i = 0; i < input_node_id_list_tmp.size(); i++)
					is_removed[input_node_id_list_tmp[i]] = false;
				// Determine the outputs
				for (int i = 0 ; i < original_module->NodeCount(); i++) {
					if (!is_removed[i]) {
						int out_edge_num = original_module->OutEdgeCount(i);
						for	(int j = 0; j < out_edge_num; j++) {
							if (is_removed[original_module->GetOutEdge(i,j)]) {
								output_node_id_list_tmp.push_back(i);
								break;
							}
						}
					}
				}
				// Add the global outputs
				for (int i = 0; i < output_node_id_list.size(); i++)
					if (!is_removed[output_node_id_list[i]]) {
						int j;
						for (j = 0; j < output_node_id_list_tmp.size(); j++)
							if (output_node_id_list[i] == output_node_id_list_tmp[j])
								break;
						if (j == output_node_id_list_tmp.size())
							output_node_id_list_tmp.push_back(output_node_id_list[i]);
					}
				// Update the input & output nodes list since the indexed will be changed after some irrelevant nodes are deleted
				for (int i = 0; i < input_node_id_list_tmp.size(); i++) {
					int c = 0;
					for (int j = 0; j < input_node_id_list_tmp[i]; j++)
						if (is_removed[j])
							c++;
					input_node_id_list_tmp[i] -= c;
				}
				for (int i = 0; i < output_node_id_list_tmp.size(); i++) {
					int c = 0;
					for (int j = 0; j < output_node_id_list_tmp[i]; j++)
						if (is_removed[j])
							c++;
					output_node_id_list_tmp[i] -= c;
				}
				// Remove irrelevant nodes
				for (int i = original_module->NodeCount() - 1; i >= 0; i--)
					if (is_removed[i])
						tmp_editor->DeleteNode(i);
				// Modify the mRNA input node to UNKNOWN
				for (int i = 0; i < input_node_id_list_tmp.size(); i++) {
					GeneCircuitComponent* tmp = ((BioNetNode*)tmp_editor->GetNodeAttr(input_node_id_list_tmp[i]))->component;
					if (tmp->getType() == m_RNA)
						tmp->name = "UNKNOWN";
				}
				Module* module_tmp = new Module(tmp_editor, input_node_id_list_tmp, output_node_id_list_tmp);
				module_tmp->name = module_name;
				//module_tmp->reference
				//cout << "=====================================" << endl;
				//module_tmp->TestPrint(&std::cout);
				//cout << "=====================================" << endl;
				module_list->push_back(module_tmp);
			}
		delete original_module;
		delete [] node_id_map;
	}
	infile.close();
	//cout << module_list->size() << endl;
	return module_list;
}

StringMatrix* PartDatabase::getAllInteractionSpecies(string species_type) {
	StringMatrix* tmp = new StringMatrix;
	QueryResultInfo qri = Init();
	qri.query = "SELECT species_name FROM species";
	qri.query += " WHERE species_type = '" + species_type + "'";
	//cout << qri.query << endl;
	StringMatrix species_name_list;
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &species_name_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	for (int i = 0; i < species_name_list.size(); i++) {
		string species_name = species_name_list[i][0];
		//cout << endl << "--------------" << endl;
		//cout << species_name << endl;
		//cout << "--------------" << endl;
		Reset(qri);
		qri.query = "SELECT src_species FROM species_interaction";
		qri.query +=  " WHERE dest_species = '" + species_name + "'";
		StringMatrix src_name_list;
		qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &src_name_list, &qri.zErrMsg);
		if(qri.rc != SQLITE_OK) {
			cout << qri.query << endl;
			std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
			sqlite3_free(qri.zErrMsg);
		}
		if (!src_name_list.empty()) {
			vector<string> list_tmp;
			list_tmp.push_back(species_name);
			for (int j = 0; j < src_name_list.size(); j++) {
				list_tmp.push_back(src_name_list[j][0]);
				//cout << src_name_list[j][0] << "\t";
			}
			tmp->push_back(list_tmp);
		}
	}
	return tmp;
}

string PartDatabase::getSpeciesType(string species_name) {
	QueryResultInfo qri = Init();
	qri.query = "SELECT species_type FROM species";
	qri.query += " WHERE species_name = '" + species_name + "'";
	//cout << qri.query << endl;
	StringMatrix regulation_type;
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &regulation_type, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (regulation_type.empty()) {
		cout << species_name << "\t is not found in the DB" <<  endl;
		return "";
	}
	else
		return regulation_type[0][0];
}
string PartDatabase::getSpeciesInteractionType(string src_species_name, string dest_species_name) {
	QueryResultInfo qri = Init();
	qri.query = "SELECT interaction_type FROM species_interaction";
	qri.query += " WHERE src_species = '" + src_species_name + "' AND dest_species = '" + dest_species_name + "'";
	//cout << qri.query << endl;
	StringMatrix interaction_type;
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &interaction_type, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (interaction_type.empty())
		return "NONE";
	else {
		//cout << src_species_name << "\t" << interaction_type[0][0] << "\t" << dest_species_name << endl;
		return interaction_type[0][0];
	}
}
int PartDatabase::getRegulationType(Molecule* m1, Molecule* m2) {
	QueryResultInfo qri = Init();
	qri.query = "SELECT regulation FROM Regulation";
	qri.query += " WHERE src_type = '" + Id2Str(m1->getType()) + "' AND src_name = '" + m1->getName() + "' AND dest_type = '" + Id2Str(m2->getType()) + "' AND dest_name = '" + m2->getName() + "'";
	//cout << qri.query << endl;
	StringMatrix regulation_type;
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &regulation_type, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	if (regulation_type.empty())
		return NONE;
	else
		if (regulation_type[0][0].compare("ACT") == STR_EQ)
			return ACTIVATORY;
		else
			return INHIBITORY;
}

int PartDatabase::getNumberOfSubsitution(string functional_module_type) const {
	QueryResultInfo qri = Init();
	qri.query = "SELECT COUNT (DISTINCT derive_id) FROM Functional_Module_Node";
	qri.query += " WHERE functional_module_type = '" + functional_module_type + "'";
	int count_tmp;
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_count, &count_tmp, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	return count_tmp;
}

Module* PartDatabase::loadSubstitution(string functional_module_type, int derive_id) const {
	QueryResultInfo qri = Init();
	ARGEdit* module_editor = new ARGEdit;
	// Read nodes
	StringMatrix module_node_list;
	Reset(qri);
	qri.query = "SELECT node_id,node_type,node_name,node_role FROM Functional_Module_Node";
	qri.query += " WHERE functional_module_type = '" + functional_module_type + "' AND derive_id = '" + Num2Str(derive_id) + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_node_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	std::map<string, int> node_id_to_bionode_id;
	IdList input_tmp, output_tmp;
	for (unsigned int node_id = 0; node_id < module_node_list.size(); node_id++) {
		int id_tmp = module_editor->InsertNode(CreateBioNetNode(module_node_list[node_id][1], module_node_list[node_id][2]));
		node_id_to_bionode_id[module_node_list[node_id][0]] = id_tmp;
		if (module_node_list[node_id][3].compare("Input") == STR_EQ)
			input_tmp.push_back(id_tmp);
		else if (module_node_list[node_id][3].compare("Output") == STR_EQ)
			output_tmp.push_back(id_tmp);
		else if (module_node_list[node_id][3].compare("Input_Output") == STR_EQ) {
			input_tmp.push_back(id_tmp);
			output_tmp.push_back(id_tmp);
		}
	}
	// Read edges
	StringMatrix module_edge_list;
	Reset(qri);
	qri.query = "SELECT src_id, dest_id, edge_type FROM Functional_Module_Edge";
	qri.query += " WHERE functional_module_type = '" + functional_module_type + "' AND derive_id = '" + Num2Str(derive_id) + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_edge_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	for (unsigned int edge_id = 0; edge_id < module_edge_list.size(); edge_id++)
		module_editor->InsertEdge(node_id_to_bionode_id[module_edge_list[edge_id][0]],node_id_to_bionode_id[module_edge_list[edge_id][1]], new BioNetEdge(module_edge_list[edge_id][2]));
	// Create a new module
	Module* tmp = new Module(module_editor);
	tmp->inputs = input_tmp;
	tmp->outputs = output_tmp;
	// TEST HERE
	//tmp->TestPrint(&std::cout);
	return tmp;
}

Module* PartDatabase::loadContent(string macro_module_type, string macro_module_name, vector<NodeVisualInfo> *module_visual_info_list) const {
	QueryResultInfo qri = Init();
	ARGEdit* module_editor = new ARGEdit;
	// Read node name
	StringMatrix node_name_list;
	Reset(qri);
	qri.query = "SELECT node_id, node_name FROM module_node";
	qri.query += " WHERE module_name = '" + macro_module_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &node_name_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	map<string, string> node_id_to_node_name;
	for (unsigned int node_id = 0; node_id < node_name_list.size(); node_id++)
		node_id_to_node_name[node_name_list[node_id][0]] = node_name_list[node_id][1];
	// Read motif name
	Reset(qri);
	qri.query = "SELECT motif_name FROM module_info";
	qri.query += " WHERE module_name = '" + macro_module_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &qri.string_matrix, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	string motif_name = qri.string_matrix[0][0];
	// Read nodes
	StringMatrix module_node_list;
	Reset(qri);
	qri.query = "SELECT node_id, node_type, node_role, pos_x, pos_y, width, height FROM motif_node";
	qri.query += " WHERE motif_name = '" + motif_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_node_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	std::map<string, int> node_id_to_bionode_id;
	IdList input_tmp, output_tmp;
	for (unsigned int node_id = 0; node_id < module_node_list.size(); node_id++) {
		if (module_visual_info_list != NULL) {
			NodeVisualInfo tmp;
			tmp.pos_x = atoi(module_node_list[node_id][3].c_str());
			tmp.pos_y = atoi(module_node_list[node_id][4].c_str());
			tmp.width = atoi(module_node_list[node_id][5].c_str());
			tmp.height = atoi(module_node_list[node_id][6].c_str());
			module_visual_info_list->push_back(tmp);
		}
		int id_tmp = module_editor->InsertNode(CreateBioNetNode(module_node_list[node_id][1], node_id_to_node_name[module_node_list[node_id][0]]));
		node_id_to_bionode_id[module_node_list[node_id][0]] = id_tmp;
		if (module_node_list[node_id][2].compare("Input") == STR_EQ)
			input_tmp.push_back(id_tmp);
		else if (module_node_list[node_id][2].compare("Output") == STR_EQ)
			output_tmp.push_back(id_tmp);
		else if (module_node_list[node_id][2].compare("Input-Output") == STR_EQ) {
			input_tmp.push_back(id_tmp);
			output_tmp.push_back(id_tmp);
		}
	}
	// Read edges
	StringMatrix module_edge_list;
	Reset(qri);
	qri.query = "SELECT src_node_id, dest_node_id, edge_type FROM motif_edge";
	qri.query += " WHERE motif_name = '" + motif_name + "'";
	qri.rc = sqlite3_exec(db, qri.query.c_str(), callback_string, &module_edge_list, &qri.zErrMsg);
	if(qri.rc != SQLITE_OK) {
		cout << qri.query << endl;
		std::cout << "SQL error: \n" << qri.zErrMsg << std::endl;
		sqlite3_free(qri.zErrMsg);
	}
	for (unsigned int edge_id = 0; edge_id < module_edge_list.size(); edge_id++)
		module_editor->InsertEdge(node_id_to_bionode_id[module_edge_list[edge_id][0]],node_id_to_bionode_id[module_edge_list[edge_id][1]], new BioNetEdge(module_edge_list[edge_id][2]));
	// Create a new module
	Module* tmp = new Module(module_editor);
	tmp->inputs = input_tmp;
	tmp->outputs = output_tmp;
	// TODO: solve for the EMPTY OR GATE here
	if (macro_module_type.compare("OR_GATE") == STR_EQ && tmp->inputs.size() == 1)
		tmp->inputs.push_back(tmp->inputs[0]);
	// TEST HERE
	//cout << macro_module_type << "\t" << macro_module_name << endl;
	//tmp->TestPrint(&std::cout);
	return tmp;
}

double ev(int k, int n, double alpha) {
	return pow(alpha,k) - 1 - n*(alpha - 1);
}
double find_scale(int k, int n) { // solve the equation (alpha^k - 1)/(alpha - 1) = n --> alpha^k - 1 - n(alpha - 1) = 0
	double u = n, l = 1 + 1e-10;
	while (u - l > 1e-10) {
		double mid = (u + l)/2;
		if ((ev(k,n,u) > 0 && ev(k,n,mid) < 0) || (ev(k,n,u) < 0 && ev(k,n,mid) > 0))
			l = mid;
		else
			u = mid;
	}
	return (u+l)/2;
}

vector<Module*>* PartDatabase::RandomizeModuleLibrary(int library_size, int max_module_size) {
	vector<Module*>* tmp;
	int total_number_of_modules = 0;
	int max_module_library_size = floor(0.4*library_size);
	double scale = find_scale(max_module_size - 1, max_module_library_size);
	for (int module_size = max_module_size; module_size > 1; module_size--) {
		int number_of_modules = round(pow(scale, max_module_size - module_size));
		total_number_of_modules += number_of_modules;
		for (int module_id = 0; module_id < number_of_modules; module_id++) {
			// Generate module topology, each module has more than one gates
			Module* module = GenerateGateNetwork(module_size,round(0.4*module_size), false);
			// Fill names
			tmp->push_back(module);
		}
	}
	// Generate single gates
	// Generate molecules cascade
	// Generate single molecules and YES-GATE
	return tmp;
}

vector<string>* split_string(string str, char c) {
	vector<string>* string_list = new vector<string>;
	int pre_index = 0;
	for (int i = 0; i < str.length(); i++)
		if (str[i] == c) {
			string_list->push_back(str.substr(pre_index, i - pre_index));
			pre_index = i + 1;
		}
	string_list->push_back(str.substr(pre_index, str.length() - pre_index));
	return string_list;
}

string ReadSentence(ifstream* f) {
	string final_sentence;
	*f >> final_sentence;
	if (final_sentence[0] == '<') {
		string tmp;
		do {
			*f >> tmp;
			final_sentence += " ";
			final_sentence += tmp;
		}
		while (tmp[tmp.length() - 1] != '>');
		final_sentence = final_sentence.erase(0,1);
		final_sentence = final_sentence.erase(final_sentence.length() - 1,1);
	}
	return final_sentence;
}

QueryResultInfo Init() {
	QueryResultInfo qri;
	qri.zErrMsg = NULL;
	qri.query = "";
	qri.string_matrix.clear();
	return qri;
}

void AppendStringMatrix(StringMatrix* sm1, StringMatrix* sm2) { // sm1 = sm1 + sm2
	for (int i = 0; i < sm2->size(); i++)
		sm1->push_back(sm2->at(i));
	delete sm2;
}

void Reset(QueryResultInfo& qri) {
	qri.zErrMsg = NULL;
	qri.query = "";
	qri.string_matrix.clear();
}

string Num2StrPlus1(int num) {
	stringstream ss;
	ss << num + 1; // PLUS 1
	string tmp(ss.str());
	return tmp;
}

string Num2Str(int num) {
	stringstream ss;
	ss << num;
	string tmp(ss.str());
	return tmp;
}
string Num2Str(double num) {
	num = round(num*1000)/1000;
	std::ostringstream strs;
	strs << num;
	return strs.str();
}
