/*
 * Framework.cpp
 *
 *  Created on: Aug 17, 2012
 *      Author: linh
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <queue>
#include <stack>
#include <iostream>
#include "Framework.h"

void OutputNodeRecord::Print() {
	//cout << "Cost: " << cost << endl;
	//cout << "Output: " << output << endl;
	PrintStringTree(&signal_list);
	//for (int i = 0; i < branch_node_info.size(); i++)
	//	cout << branch_node_info[i].signal_name << ((branch_node_info[i].is_contained)? " Y ":" N ");
	//cout << endl;
	//cout << "Module id: " << module_id << endl;
	//cout << "Matching id: " << matching_id << endl;
	//cout << "Previous nodes: ";
	//for (int i = 0; i < pre_node_record_list.size(); i++)
	//	cout << pre_node_record_list[i] << " ";
	//cout << endl;
}

Framework::Framework(PartDatabase* partDB_): partDB(partDB_) {
	// Do nothing
}

Framework::Framework(PartDatabase* partDB_, DesignProblem dp_): partDB(partDB_), dp(dp_){
	// Do nothing
}

Framework::~Framework() {
	// TODO Auto-generated destructor stub
}

vector<GeneCircuitGraph*>* Framework::ModuleMatch(int number_of_solutions, int synthetic_module_library_size, int number_of_DB_replicates) {
	vector<GeneCircuitGraph*>* solution_list = new vector<GeneCircuitGraph*>; // TODO: clean
	vector<Module*> *module_lib; // TODO: clean
	if (synthetic_module_library_size > 1)
		module_lib = partDB->RandomizeModuleLibrary(synthetic_module_library_size, 15);
	else
		module_lib = partDB->loadModuleLibrary();
	// TEST HERE
	//dp.gcg.TestPrint(&std::cout);
	// END TEST

	ARGEdit* aux_editor = (ARGEdit*) dp.gcg.ConvertToLoader(); // TODO: clean
	// Insert linkers
	int n = dp.gcg.NodeCount();
	IdList linker_node_id_list;
	for (int node_id = 0; node_id < n; node_id++) {
		dp.gcg.GetNodeAttr(node_id)->original_node_id = node_id;
		((BioNetNode*)aux_editor->GetNodeAttr(node_id))->original_node_id = node_id;
		if (dp.gcg.GetNodeAttr(node_id)->component->getCategory() == LOGIC_GATE || (dp.gcg.GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL && dp.gcg.GetNodeAttr(node_id)->component->getType() != LINKER) ) {
			int adj_n = dp.gcg.OutEdgeCount(node_id);
			bool is_connected_to_output = false;
			for (int adj_i = 0; adj_i < adj_n; adj_i++) {
				GeneCircuitComponent *comp = dp.gcg.GetNodeAttr(dp.gcg.GetOutEdge(node_id, adj_i))->component;
				if (comp->getCategory() == MOLECULE) {
					is_connected_to_output = true;
					break;
				}
			}
			if (!is_connected_to_output && adj_n > 0)
				linker_node_id_list.push_back(node_id);
		}
	}
	for (int i = 0; i < linker_node_id_list.size(); i++) {
		int node_id = linker_node_id_list[i];
		int adj_n = dp.gcg.OutEdgeCount(node_id);
		int id_tmp = aux_editor->InsertNode(CreateBioNetNode("LINKER"));
		//((BioNetNode*)aux_editor->GetNodeAttr(id_tmp))->original_node_id = node_id;
		aux_editor->InsertEdge(node_id, id_tmp, new BioNetEdge("UNKNOWN"));
		for (int adj_i = 0; adj_i < adj_n; adj_i++) {
			int adj_v_i = dp.gcg.GetOutEdge(node_id, adj_i);
			aux_editor->InsertEdge(id_tmp, adj_v_i, dp.gcg.GetEdgeAttr(node_id, adj_v_i)->clone());
			aux_editor->DeleteEdge(node_id, adj_v_i);
		}
	}
	GeneCircuitGraph* extension_gcg = new GeneCircuitGraph(aux_editor); // TODO: clean
	// TEST HERE
	extension_gcg->TestPrint(&std::cout);
	// END TEST

	// Match each module with the input graph to build the matching map
	GraphIsomorphismCollection** matching_map = new GraphIsomorphismCollection*[module_lib->size()]; // TODO: clean
	vector<IdPairList> matching_list_at_node;
	matching_list_at_node.resize(extension_gcg->NodeCount());
	for (unsigned int module_id = 0; module_id < module_lib->size(); module_id++) {
		// call Ullman here
		Module* module_tmp = module_lib->at(module_id); // TODO: clean
		matching_map[module_id] = new GraphIsomorphismCollection;
		VF2SubState vf2_sub_state(module_tmp, extension_gcg);
		match(&vf2_sub_state, all_visitor, matching_map[module_id]); // TODO: we should sort pairs in each isomorphism so that it is indexed by the small node id
		// TEST HERE
		//cout << "---------------" << endl;
		//module_tmp->TestPrint(&std::cout);
		//cout << "---------------" << endl;
		//for (unsigned int matching_id = 0; matching_id < matching_map[module_id]->size(); matching_id++) {
		//	cout << "Matching: " << matching_id + 1 << "-th" << endl;
		//	for (unsigned int node_pair_id = 0; node_pair_id < matching_map[module_id]->at(matching_id).size(); node_pair_id++)
		//		cout << matching_map[module_id]->at(matching_id)[node_pair_id].large_graph_node << "\t" << matching_map[module_id]->at(matching_id)[node_pair_id].small_graph_node << endl;
		//}
		// END TEST
		bool *is_input = new bool[module_tmp->NodeCount()];
		for (int i = 0; i < module_tmp->NodeCount(); i++)
			is_input[i] = false;
		for (int i = 0; i < module_tmp->inputs.size(); i++)
			is_input[module_tmp->inputs[i]] = true;
		for (unsigned int matching_id = 0; matching_id < matching_map[module_id]->size(); matching_id++) {
			if (module_tmp->outputs.empty()) {
				cout << "SBROME: Error, module has no output" << endl;
				module_tmp->TestPrint(&std::cout);
			}
			int module_node_output_index = module_tmp->outputs[0]; // TODO: extend for the case of multiple outputs
			// NOTICE:
			//		+ matching_list_at_node[node_id][i].first: id of the matched module_id
			//		+ matching_list_at_node[node_id][i].second: id of the matching corresponded to module_id

			bool* is_matched = new bool[extension_gcg->NodeCount()];	// use to check if a node is matched with another node in the module graph
			for (int node_id = 0; node_id < extension_gcg->NodeCount(); node_id++)
				is_matched[node_id] = false;
			for (unsigned int node_pair_id = 0; node_pair_id < matching_map[module_id]->at(matching_id).size(); node_pair_id++)
				is_matched[matching_map[module_id]->at(matching_id)[node_pair_id].large_graph_node] = true;
			bool encapsulated_checking = true;	// check to ensure that there is no connection between a node outside of the module with a non-input node inside the module
			for (unsigned int node_pair_id = 0; node_pair_id < matching_map[module_id]->at(matching_id).size(); node_pair_id++) {
				int circuit_node_id = matching_map[module_id]->at(matching_id)[node_pair_id].large_graph_node;
				int in_degree = extension_gcg->InEdgeCount(circuit_node_id);
				for (int in_node_id = 0; in_node_id < in_degree; in_node_id++)
					if (!is_matched[extension_gcg->GetInEdge(circuit_node_id,in_node_id)] && !is_input[matching_map[module_id]->at(matching_id)[node_pair_id].small_graph_node]) {
						encapsulated_checking = false;
						break;
					}
				if (!encapsulated_checking)
					break;
			}
			if (encapsulated_checking)
				matching_list_at_node[matching_map[module_id]->at(matching_id)[module_node_output_index].large_graph_node].push_back(std::make_pair(module_id,matching_id));
			delete []is_matched;
		}
		delete []is_input;
	}

	// Labeling nodes
	int number_of_nodes = extension_gcg->NodeCount();
	bool* marked = new bool[number_of_nodes];	// TODO: clean
	int* in_deg = new int[number_of_nodes];		// TODO: clean // in-degree list of nodes
	int* topo_order = new int[number_of_nodes];	// TODO: clean
	for (int i = 0; i < number_of_nodes; i++) {
		in_deg[i] = extension_gcg->InEdgeCount(i);
		marked[i] = false;
	}
	int current_index = 0;
	do {
		for (int i = 0; i < number_of_nodes; i++)
			if (!marked[i] && in_deg[i] == 0) {
				int out_num = extension_gcg->OutEdgeCount(i);
				for (int j = 0; j < out_num; j++)
					in_deg[extension_gcg->GetOutEdge(i,j)]--;
				marked[i] = true;
				topo_order[current_index] = i;
				current_index++;
			}
	}
	while (current_index < number_of_nodes);
	delete [] marked;
	delete [] in_deg;

	// Determine the branching nodes
	map<int,int> branching_node_map;
	int number_of_branching_nodes = 0;
	for (int node_id = 0; node_id < number_of_nodes; node_id++)
		if (extension_gcg->OutEdgeCount(node_id) > 1)
			branching_node_map[node_id] = number_of_branching_nodes++;

	// Main loop here to update the information record for each intermediate output node
	vector<OutputNodeRecordList> record_list;
	record_list.resize(number_of_nodes);
	for (int topo_id = 0; topo_id < number_of_nodes; topo_id++) {
		// TEST HERE
		// cout << topo_order[i] << endl;
		// END TEST
		int node_id = topo_order[topo_id];
		if (extension_gcg->GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL || extension_gcg->GetNodeAttr(node_id)->component->getCategory() == MOLECULE) { // TODO: extend to "select" node for the future
			if (extension_gcg->GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL && extension_gcg->GetNodeAttr(node_id)->component->getType() == INPUT) {	// input nodes
				record_list[node_id].resize(1);
				record_list[node_id][0].output = "UNKNOWN";
				record_list[node_id][0].branch_node_info.resize(number_of_branching_nodes);
			}
			else {
				for (unsigned int matching_id = 0; matching_id < matching_list_at_node[node_id].size(); matching_id++) {
					int matched_module_id = matching_list_at_node[node_id][matching_id].first;
					Module* current_matched_module = module_lib->at(matched_module_id);
					int matching_id_in_mapping = matching_list_at_node[node_id][matching_id].second;
					// TEST HERE
					//cout << "=================" << endl;
					//cout << "Node: " << node_id << "\t Module:" << module_lib->at(matched_module_id)->name << endl;
					//cout << "=================" << endl;
					//for (unsigned int i = 0; i < matching_map[matched_module_id]->at(matching_id_in_mapping).size(); i++) {
					//	cout << matching_map[matched_module_id]->at(matching_id_in_mapping)[i].large_graph_node << "\t" << matching_map[matched_module_id]->at(matching_id_in_mapping)[i].small_graph_node << endl;				}
					// END TEST
					// main loop to compute the product of all fan-ins
					int number_of_fan_ins = current_matched_module->inputs.size();
					int* fan_in_node_id_list = new int[number_of_fan_ins];
					// TODO: clear
					int* input_node_id_list = new int[number_of_fan_ins];
					// TODO: clear
					int* numbering_array = new int[number_of_fan_ins];
					// TODO: clear
					bool record_is_empty = false;
					for (int fan_index = 0; fan_index < number_of_fan_ins; fan_index++) {
						fan_in_node_id_list[fan_index] = matching_map[matched_module_id]->at(matching_id_in_mapping)[current_matched_module->inputs[fan_index]].large_graph_node;
						input_node_id_list[fan_index] = current_matched_module->inputs[fan_index];
						// TEST HERE
						//cout << fan_in_node_id_list[fan_index] << endl;
						// END TEST
						if (record_list[fan_in_node_id_list[fan_index]].size() > 0)
							numbering_array[fan_index] = 0;
						else {
							record_is_empty = true;
							break;
						}
					}
					if (!record_is_empty) {
						int idx_tmp;
						do {
							// Synthesize a new record here
							OutputNodeRecord composite_record_tmp;
							composite_record_tmp.branch_node_info.resize(number_of_branching_nodes);
							composite_record_tmp.cost = 1;
							StringList empty_str_list;
							bool is_output_compatible = true;
							bool is_cross_talking = false;
							bool is_branching_node_ok = true;
							// Update the signal check list and the branching node check list with the mapping at the output node
							for(int node_id_tmp = 0; node_id_tmp < current_matched_module->NodeCount(); node_id_tmp++){
								int node_index_in_input_graph = matching_map[matched_module_id]->at(matching_id_in_mapping)[node_id_tmp].large_graph_node;
								if (branching_node_map.count(node_index_in_input_graph) > 0) {
									composite_record_tmp.branch_node_info[branching_node_map[node_index_in_input_graph]].signal_name = current_matched_module->GetNodeAttr(node_id_tmp)->component->name;
									composite_record_tmp.branch_node_info[branching_node_map[node_index_in_input_graph]].is_contained = (current_matched_module->InEdgeCount(node_id_tmp) > 1 || extension_gcg->InEdgeCount(node_index_in_input_graph) == 0);
								}
							}
							// Update these two check lists with the mapping at each fan-in
							IdList temp_update_node_list;	// contains nodes in the matched module where name is UNKNOWN, it is updated temporarily with the name of the fan-in module output
							for (int fan_index = 0; fan_index < number_of_fan_ins; fan_index++) {
								OutputNodeRecord record_tmp = record_list[fan_in_node_id_list[fan_index]][numbering_array[fan_index]];
								composite_record_tmp.cost += record_tmp.cost;
								// Check compatible
								is_output_compatible = (record_tmp.output.compare("UNKNOWN") == STR_EQ || current_matched_module->GetNodeAttr(input_node_id_list[fan_index])->component->name.compare("UNKNOWN") == STR_EQ || current_matched_module->GetNodeAttr(input_node_id_list[fan_index])->component->name.compare(record_tmp.output) == STR_EQ);
								if (!is_output_compatible) {
									//cout << "NOT COMPATIBLE, output: " << record_tmp.output << "\t input needed: " << current_matched_module->GetNodeAttr(input_node_id_list[fan_index])->component->name << endl;
									break;
								}
								else if (current_matched_module->GetNodeAttr(input_node_id_list[fan_index])->component->name.compare("UNKNOWN") == STR_EQ) {
									temp_update_node_list.push_back(input_node_id_list[fan_index]);
									current_matched_module->GetNodeAttr(input_node_id_list[fan_index])->component->name = record_tmp.output;
								}
								// Check cross-talk
								// if the fan-in signal output is connected to a pool node, skip the cross-talk check on this signal
								int pool_node_id = current_matched_module->GetOutEdge(input_node_id_list[fan_index],0);
								vector<Molecule*>* mol_list_tmp = SrcComponent2Molecule(current_matched_module->GetNodeAttr(pool_node_id)->component, partDB);
								CheckList skip_list;
								if (mol_list_tmp->at(0)->getCategory() == MOLECULE && mol_list_tmp->at(0)->getType() == POOL)
									skip_list[record_tmp.output] = 1;
								for (int i = 0; i < mol_list_tmp->size(); i++)
									delete mol_list_tmp->at(i);
								delete mol_list_tmp;
								is_cross_talking = check_and_merge(&composite_record_tmp.signal_list, &record_tmp.signal_list, &skip_list);
								if (is_cross_talking) {
									//cout << "CROSS TALKING" << endl;
									break;
								}
								// Check branching-nodes
								for (int branching_node_id = 0; branching_node_id < number_of_branching_nodes; branching_node_id++) {
									if (!record_tmp.branch_node_info[branching_node_id].signal_name.empty())
										if (composite_record_tmp.branch_node_info[branching_node_id].signal_name.empty() || (composite_record_tmp.branch_node_info[branching_node_id].signal_name.compare(record_tmp.branch_node_info[branching_node_id].signal_name) == STR_EQ))
											composite_record_tmp.branch_node_info[branching_node_id] = record_tmp.branch_node_info[branching_node_id];
										else {
											is_branching_node_ok = false;
											break;
										}
								}
								if (!is_branching_node_ok) {
									//cout << "BRANCHING NODE IS NOT OK" << endl;
									break;
								}
							}
							if (is_output_compatible && !is_cross_talking && is_branching_node_ok) {
								StrTree matched_module_str_tree;
								GCG2StrTree(Layout(current_matched_module, partDB), &matched_module_str_tree);
								// TEST HERE
								//cout << "::::::::::::::" << endl;
								//PrintStringTree(&matched_module_str_tree);
								//cout << ";;;;;;;;;;;;;;" << endl;
								//PrintStringTree(&composite_record_tmp.signal_list);
								//cout << "::::::::::::::" << endl;
								// END TEST
								if (!check_and_merge(&composite_record_tmp.signal_list, &matched_module_str_tree, NULL)) {
									composite_record_tmp.output = current_matched_module->getOutputName()[0];
									composite_record_tmp.module_id = matched_module_id;
									composite_record_tmp.matching_id = matching_id_in_mapping;
									composite_record_tmp.pre_node_record_list.resize(number_of_fan_ins);
									for (int i = 0; i < number_of_fan_ins; i++)
										composite_record_tmp.pre_node_record_list[i] = numbering_array[i];
									// Add this new record
									record_list[node_id].push_back(composite_record_tmp);
									// TEST HERE
									//cout << "^-^" << endl;
									// END TEST
								}
								//else
									//cout << "CROSS TALKING" << endl;
							}
							// Recover some input nodes of the module that we updated temporarily back to UNKNOWN
							for (int id_tmp = 0; id_tmp < temp_update_node_list.size(); id_tmp++)
								current_matched_module->GetNodeAttr(temp_update_node_list[id_tmp])->component->name = "UNKNOWN";
							// Re-numbering
							idx_tmp = 0;
							while ((idx_tmp < number_of_fan_ins) && (numbering_array[idx_tmp] == record_list[fan_in_node_id_list[idx_tmp]].size() - 1)) {
								numbering_array[idx_tmp] = 0;
								idx_tmp++;
							}
							if (idx_tmp < number_of_fan_ins)
								numbering_array[idx_tmp]++;
						}
						while (idx_tmp < number_of_fan_ins);
					}
				}
			}
		}
		// TEST HERE
		//cout << "****************************" << endl;
		//cout << "Node: " << node_id << endl;
		//cout << "****************************" << endl;
		//for (int i = 0; i < record_list[node_id].size(); i++) {
		//	cout << "+++++++++++++++++++++++++ record " << i << "-th" << endl;
		//	record_list[node_id][i].Print();
		//}
		// END TEST
	}

	// Trace back to find the optimal solution
	int output_node_id = topo_order[number_of_nodes - 1];
	int optimal_index = 0;
	if (record_list[output_node_id].size() > 0) {
		int optimal_cost = record_list[output_node_id][0].cost;
		for (int i = 1; i < record_list[output_node_id].size(); i++)
			if (record_list[output_node_id][i].cost < optimal_cost) {
				optimal_cost = record_list[output_node_id][i].cost;
				optimal_index = i;
			}
		queue<IdPair> traverse_queue;	// each pair contains the output node id (in the large graph) and the record id (i.e. the optimal record for this node)
		traverse_queue.push(std::make_pair(output_node_id, optimal_index));
		GeneCircuitGraph* result_gcg = extension_gcg->clone();
		while (!traverse_queue.empty()){
			IdPair pair_tmp = traverse_queue.front();
			traverse_queue.pop();
			int matched_module_id = record_list[pair_tmp.first][pair_tmp.second].module_id;
			int matching_id = record_list[pair_tmp.first][pair_tmp.second].matching_id;
			Module* matched_module = module_lib->at(matched_module_id);
			// TEST HERE
			//cout << "Node id: " << pair_tmp.first << "\tRecord id: " << pair_tmp.second << "\tModule id: " << matched_module_id << "\tMatching id: " << matching_id << "\tTotal matches: " << matching_map[matched_module_id]->size() << endl;
			// END TEST
			// Assign the name to the result graph nodes
			for (int i = 0; i < matching_map[matched_module_id]->at(matching_id).size(); i++) {
				int result_graph_node_id = matching_map[matched_module_id]->at(matching_id)[i].large_graph_node;
				int module_graph_node_id = matching_map[matched_module_id]->at(matching_id)[i].small_graph_node;
				if (matched_module->InEdgeCount(module_graph_node_id) > 0 || result_gcg->InEdgeCount(result_graph_node_id) == 0) {	// i.e., an input node of a modules
					result_gcg->GetNodeAttr(result_graph_node_id)->component->name = matched_module->GetNodeAttr(module_graph_node_id)->component->name;
					if (result_gcg->GetNodeAttr(result_graph_node_id)->component->getCategory() == CIRCUIT_SIGNAL)
						if (matched_module->GetNodeAttr(module_graph_node_id)->component->getCategory() == CIRCUIT_SIGNAL)
							((CircuitSignal*)result_gcg->GetNodeAttr(result_graph_node_id)->component)->physical_instance = ((CircuitSignal*)matched_module->GetNodeAttr(module_graph_node_id)->component)->physical_instance->clone();
						else
							((CircuitSignal*)result_gcg->GetNodeAttr(result_graph_node_id)->component)->physical_instance = matched_module->GetNodeAttr(module_graph_node_id)->component->clone();
				}
			}
			// push the fan-ins into the queue
			int number_of_fan_ins = module_lib->at(matched_module_id)->inputs.size();
			for (int i = 0; i < number_of_fan_ins; i++) {
				int fan_in_node_id_in_result_graph = matching_map[matched_module_id]->at(matching_id)[module_lib->at(matched_module_id)->inputs[i]].large_graph_node;
				if (result_gcg->GetNodeAttr(fan_in_node_id_in_result_graph)->component->getCategory() != CIRCUIT_SIGNAL || result_gcg->GetNodeAttr(fan_in_node_id_in_result_graph)->component->getType() != INPUT) // i.e., not an input node of the result graph
					traverse_queue.push(std::make_pair(fan_in_node_id_in_result_graph,record_list[pair_tmp.first][pair_tmp.second].pre_node_record_list[i]));
			}
		}
		// TEST HERE
		result_gcg->TestPrint(&std::cout);
		// END TEST
		solution_list->push_back(result_gcg);
	}

	delete [] topo_order;
	delete extension_gcg;
	return solution_list;
}

void GCG2StrTree(GeneCircuitGraph* gcg, StrTree* str_tree) {
	int n = gcg->NodeCount();
	for (int node_id = 0; node_id < gcg->NodeCount(); node_id++) {
		int in_deg = gcg->InEdgeCount(node_id);
		if (gcg->GetNodeAttr(node_id)->component->getCategory() == MOLECULE && in_deg > 0) {	// TODO: fix for the pool nodes also
			StringList str_list_tmp;
			for (int i = 0; i < in_deg; i++)
				str_list_tmp.push_back(gcg->GetNodeAttr(gcg->GetInEdge(node_id,i))->component->name);
			(*str_tree)[gcg->GetNodeAttr(node_id)->component->name].first = str_list_tmp;
			(*str_tree)[gcg->GetNodeAttr(node_id)->component->name].second = -1;	// is not the global input
		}
	}
	// TEST HERE
	//for (StrTree::iterator it = str_tree->begin(); it != str_tree->end(); it++) {
	//	cout << it->first << "  ";
	//	for (int i = 0; i < it->second.size(); i++)
	//		cout << " --" << it->second[i];
	//	cout << endl;
	//}
	// END TEST
}

bool equivalence_check(StrTree* original_str_tree, StrTree* sub_str_tree, string signal) {
	// TEST HERE
	//cout << "************" << endl;
	//cout << signal << endl;
	//PrintStringTree(original_str_tree);
	//PrintStringTree(sub_str_tree);
	//cout << "************" << endl;
	// END TEST HERE
	StringList sl1;
	if (original_str_tree->count(signal))
		sl1 = original_str_tree->at(signal).first;
	StringList sl2;
	if (sub_str_tree->count(signal))
		sl2 = sub_str_tree->at(signal).first;
	if (sl1.size() != sl2.size())
		return false;
	else
		for (int i = 0; i < sl1.size(); i++) {
			int j;
			for (j = 0; j < sl2.size(); j++)
				if (sl1[i].compare(sl2[j]) == STR_EQ) {
					if (!equivalence_check(original_str_tree, sub_str_tree, sl1[i]))
						return false;
					break;
				}
			if (j == sl2.size())
				return false;
		}
	return true;
}

bool check_and_merge(StrTree* original_str_tree, StrTree* sub_str_tree, CheckList* skip_list) {
	for (StrTree::iterator it = sub_str_tree->begin(); it != sub_str_tree->end(); ++it)
		if (original_str_tree->count(it->first) == 0)
			(*original_str_tree)[it->first] = it->second;
		else {
			if (it->second.second >= 0 && (*original_str_tree)[it->first].second >= 0 && it->second.second != (*original_str_tree)[it->first].second) {
				//cout << "Cross-talk with: " << it->first << endl;
				return true;
			}
			if (!equivalence_check(original_str_tree, sub_str_tree, it->first))
				if(skip_list == NULL || skip_list->count(it->first) == 0) {
					//cout << "Cross-talk with: " << it->first << endl;
					return true;
				}
				else
					UnionStrList(&(*original_str_tree)[it->first].first, &it->second.first);
		}
	return false;	// it means there is no cross-talking
}

void UnionStrList(StringList* first_str_list, StringList* second_str_list) {
	int l1 = first_str_list->size();
	int l2 = second_str_list->size();
	for (int j = 0; j < l2; j++) {
		int i;
		for (i = 0; i < l1; i++)
			if (first_str_list->at(i).compare(second_str_list->at(j)) == STR_EQ)
				break;
		if (i == l1)
			first_str_list->push_back(second_str_list->at(j));
	}
}

void PrintStringTree(StrTree* str_tree) {
	cout << "~~~~~~~~~~" << endl;
	for (StrTree::iterator it = str_tree->begin(); it != str_tree->end(); it++) {
		cout << it->first << ((it->second.second >= 0)?"(Input) ":" ");
		for (int i = 0; i < it->second.first.size(); i++)
			cout << " --" << it->second.first[i];
		cout << endl;
	}
	cout << "~~~~~~~~~~" << endl;
}

vector<Molecule*>* SrcComponent2Molecule(GeneCircuitComponent *component, PartDatabase* partDB) {
	vector<Molecule*>* tmp = new vector<Molecule*>;
	if (component->getCategory() == MOLECULE)
		tmp->push_back((Molecule*) component->clone());
	else if (component->getCategory() == CIRCUIT_SIGNAL)
		tmp->push_back((Molecule*)((CircuitSignal*) component)->physical_instance->clone());
	else {
		Module* module_tmp = partDB->loadContent(Id2Str(component->getType()),component->name);
		// TEST HERE
		//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ " << Id2Str(component->getType()) << "\t"<< component->name << endl;
		//module_tmp->TestPrint(&std::cout);
		//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		// END TEST
		for (unsigned int i = 0; i < module_tmp->inputs.size(); i++) {
			vector<Molecule*>* molecule_list_tmp  = SrcComponent2Molecule(module_tmp->GetNodeAttr(module_tmp->inputs[i])->component, partDB);
			for (unsigned int j = 0; j < molecule_list_tmp->size(); j++)
				tmp->push_back((Molecule*)molecule_list_tmp->at(j)->clone());
			delete molecule_list_tmp;
		}
		delete module_tmp;
	}
	return tmp;
}

Molecule* DestComponent2Molecule(GeneCircuitComponent *component, PartDatabase* partDB) {
	if (component->getCategory() == MOLECULE)
		return ((Molecule*) component->clone());
	else if (component->getCategory() == CIRCUIT_SIGNAL)
		return ((Molecule*)((CircuitSignal*)component)->physical_instance->clone());
	else {
		Module* module_tmp = partDB->loadContent(Id2Str(component->getType()),component->name);
		GeneCircuitComponent* component_tmp = DestComponent2Molecule(module_tmp->GetNodeAttr(module_tmp->outputs[0])->component, partDB)->clone();
		delete module_tmp;
		return ((Molecule*)component_tmp);
	}
}

int Find_Regulation(GeneCircuitComponent* src, GeneCircuitComponent* dest, PartDatabase *partDB) {
	int regulation_type = UNKNOWN;
	if (src->name.length() > 0 && dest->name.length() > 0) {
		Molecule* src_mol = DestComponent2Molecule(src, partDB);
		//cout << "src: " << src_mol->name << endl;
		vector<Molecule*>* dest_mol_list = SrcComponent2Molecule(dest, partDB);
		for (unsigned int i = 0; i < dest_mol_list->size(); i++) {
			int t = partDB->getRegulationType(src_mol, dest_mol_list->at(i));
			//cout << "dest " << i << ":\t" <<  dest_mol_list->at(i)->name << "\t" << t << endl;
			if (t != NONE) {
				regulation_type = t;
				break;
			}
		}
		if (regulation_type == UNKNOWN)
			regulation_type = NONE;
		for (unsigned int i = 0; i < dest_mol_list->size(); i++)
			delete dest_mol_list->at(i);
		delete dest_mol_list;
	}
	return regulation_type;
}

bool Check_Update(GeneCircuitGraph *gcg, int src_node_id, int dest_node_id, PartDatabase* partDB) {
	GeneCircuitComponent *src = gcg->GetNodeAttr(src_node_id)->component;
	GeneCircuitComponent *dest = gcg->GetNodeAttr(dest_node_id)->component;
	if (src->name.length() > 0 && dest->name.length() > 0) {
		Molecule* src_mol = DestComponent2Molecule(src, partDB);
		vector<Molecule*>* dest_mol_list = SrcComponent2Molecule(dest, partDB);
		bool compatible = false;
		for (unsigned int i = 0; i < dest_mol_list->size(); i++) {
			int type_1 = partDB->getRegulationType(src_mol, dest_mol_list->at(i));
			int type_2 = gcg->GetEdgeAttr(src_node_id, dest_node_id)->type;
			// TEST HERE
			//cout << src_mol->name << " " << type_1 << " " << dest_mol_list->at(i)->name
			//	<< "\t" << src_node_id << " " << type_2 << " " << dest_node_id << endl;
			if ((type_1 == type_2) || (type_1 != NONE && type_2 == UNKNOWN)) {
				compatible = true;
				break;
			}
		}
		delete dest_mol_list;
		return compatible;
	}
	// TODO: update src_mol and dest_mol
	return true;
}

GeneCircuitGraph* Layout(GeneCircuitGraph *gcg, PartDatabase* partDB, map<int,NodeVisualInfo*>* visual_info_list) {
	// Remove linkers first
	ARGEdit* editor = (ARGEdit*) gcg->ConvertToLoader();
	IdList deleted_node_list;
	for (int node_id = 0; node_id < gcg->NodeCount(); node_id++) {
		//if (gcg->GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL && gcg->GetNodeAttr(node_id)->component->getType() == LINKER) {
		if (gcg->GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL) {
			if (visual_info_list == NULL || (gcg->GetNodeAttr(node_id)->component->getType() == LINKER))
				deleted_node_list.push_back(node_id);
			if (gcg->GetNodeAttr(node_id)->component->getType() == LINKER) {
				int src_id = gcg->GetInEdge(node_id,0);
				int number_of_succ_nodes = gcg->OutEdgeCount(node_id);
				for (int succ_id = 0; succ_id < number_of_succ_nodes; succ_id++) {
					int succ_node_id = gcg->GetOutEdge(node_id, succ_id);
					editor->InsertEdge(src_id, succ_node_id, gcg->GetEdgeAttr(node_id, succ_node_id)->clone());
				}
			}
		}
	}
	for (unsigned int i = 0; i < deleted_node_list.size(); i++)
		editor->DeleteNode(deleted_node_list[i] - i);
	// Extend device nodes
	GeneCircuitGraph *current_gcg = new GeneCircuitGraph(editor);
	bool continue_extending;
	do {
		continue_extending = false;
		for (int node_id = 0; node_id < current_gcg->NodeCount(); node_id++)
			if (current_gcg->GetNodeAttr(node_id)->component->getCategory() == CIRCUIT_SIGNAL) {
				//GeneCircuitComponent* comp_tmp;
				//if (current_gcg->GetNodeAttr(node_id)->component->getType() == OUTPUT)
				//	comp_tmp = CreateComponent("Pool", current_gcg->GetNodeAttr(node_id)->component->name);
				//else
				//	comp_tmp = ((CircuitSignal*)current_gcg->GetNodeAttr(node_id)->component)->physical_instance->clone();
				//delete current_gcg->GetNodeAttr(node_id)->component;
				//current_gcg->GetNodeAttr(node_id)->component = comp_tmp;
			}
			else if (current_gcg->GetNodeAttr(node_id)->component->getCategory() != MOLECULE && !current_gcg->GetNodeAttr(node_id)->component->getName().empty()) { // check if the name is empty since we can keep the device node by setting it's name = empty
				vector<NodeVisualInfo> *module_visual_info_list = new vector<NodeVisualInfo>;
				Module *module_tmp = partDB->loadContent(Id2Str(current_gcg->GetNodeAttr(node_id)->component->getType()), current_gcg->GetNodeAttr(node_id)->component->name, module_visual_info_list);
				GeneCircuitGraph *gcg_tmp = GeneCircuitGraphExtension(current_gcg, node_id, module_tmp, partDB, visual_info_list, module_visual_info_list);
				delete current_gcg;
				current_gcg = gcg_tmp;
				continue_extending = true;
				break;
			}
		// TEST HERE
		//cout << "***************************";
		//current_gcg->TestPrint(&std::cout);
	} while (continue_extending);
	// Update the edges that are created during the expansion step
	for (int src_id = 0; src_id < current_gcg->NodeCount(); src_id++)
		for (int dest_id = 0; dest_id < current_gcg->NodeCount(); dest_id++)
			if (src_id != dest_id) {
				BioNetEdge *bn_e = current_gcg->GetEdgeAttr(src_id, dest_id);
				if (bn_e != NULL && bn_e->type == UNKNOWN)
					bn_e->type = Find_Regulation(current_gcg->GetNodeAttr(src_id)->component, current_gcg->GetNodeAttr(dest_id)->component, partDB);
			}
	return current_gcg;
}

GeneCircuitGraph* GeneCircuitGraphExtension(GeneCircuitGraph* gcg, int sub_node_id, Module* module, PartDatabase *partDB, map<int,NodeVisualInfo*>* visual_info_list, vector<NodeVisualInfo> *module_visual_info_list) {
	// TEST HERE
	//cout << "---------------------------";
	//gcg->TestPrint(&std::cout);
	//cout << "---------";
	//module->TestPrint(&std::cout);
	//cout << "---------------------------";
	// END TEST
	ARGEdit *editor_tmp = new ARGEdit;
	map<int,int> circuit_to_extension;
	// Add nodes
	for (int i = 0; i < gcg->NodeCount(); i++)
		if (i != sub_node_id)
			circuit_to_extension[i] = editor_tmp->InsertNode(gcg->GetNodeAttr(i)->clone());
		else if (visual_info_list != NULL) { // visual_info_list != NULL ---> keep the expanded node for the visual expansion
			gcg->GetNodeAttr(i)->component->name = "";	// set the name to be empty so that it will not be expanded later
			circuit_to_extension[i] = editor_tmp->InsertNode(gcg->GetNodeAttr(i)->clone());
		}
	// Add edges
	for (int i = 0; i < gcg->NodeCount(); i++)
		for (int j = 0; j < gcg->NodeCount(); j++)
			if (i != sub_node_id && j != sub_node_id && i != j && gcg->GetEdgeAttr(i,j) != NULL) // all edges connect to the expanded node the will be removed
				editor_tmp->InsertEdge(circuit_to_extension[i], circuit_to_extension[j], gcg->GetEdgeAttr(i,j)->clone());
	map<int,int> module_to_extension;
	// Add nodes
	double x_min = 10000, y_min = 10000, x_max = 0, y_max = 0;
	double w_z, h_z; // for the zone of this motif
	if (module_visual_info_list != NULL) {
		for (unsigned int node_id = 0; node_id < module_visual_info_list->size(); node_id++) {
			if (module_visual_info_list->at(node_id).pos_x < x_min)
				x_min = module_visual_info_list->at(node_id).pos_x;
			if (module_visual_info_list->at(node_id).pos_x + module_visual_info_list->at(node_id).width > x_max)
				x_max = module_visual_info_list->at(node_id).pos_x + module_visual_info_list->at(node_id).width;
			if (module_visual_info_list->at(node_id).pos_y < y_min)
				y_min = module_visual_info_list->at(node_id).pos_y;
			if (module_visual_info_list->at(node_id).pos_y + module_visual_info_list->at(node_id).height > y_max)
				y_max = module_visual_info_list->at(node_id).pos_y + module_visual_info_list->at(node_id).height;
		}
		w_z = x_max - x_min;
		h_z = y_max - y_min;
	}
	double x_parent, y_parent, w_parent, h_parent;	// visual information of the parent node
	double x_c, y_c, w_c, h_c; 	// visual information of the container
	if (visual_info_list != NULL) {
		x_parent = visual_info_list->at(sub_node_id)->pos_x;	// visual information of the parent node
		y_parent = visual_info_list->at(sub_node_id)->pos_y;
		w_parent = visual_info_list->at(sub_node_id)->width;
		h_parent = visual_info_list->at(sub_node_id)->height;
		switch (gcg->GetNodeAttr(sub_node_id)->component->getType()) {
			case YES_GATE:
				x_c = x_parent + 0.03*w_parent;
				y_c = y_parent + 0.25*h_parent;
				w_c = 0.8*w_parent;
				h_c = 0.5*h_parent;
				break;
			case NOT_GATE:
				x_c = x_parent + 0.03*w_parent;
				y_c = y_parent + 0.25*h_parent;
				w_c = 0.8*w_parent;
				h_c = 0.5*h_parent;
				break;
			case AND_GATE:
				x_c = x_parent + 0.05*w_parent;
				y_c = y_parent + 0.1*h_parent;
				w_c = 0.9*w_parent;
				h_c = 0.8*h_parent;
				break;
			case OR_GATE:
				x_c = x_parent + 0.5*w_parent;
				y_c = y_parent + 0.25*h_parent;
				w_c = 0.5*w_parent;
				h_c = 0.5*h_parent;
				break;
		}
	}
	for (int i = 0; i < module->NodeCount(); i++) {
		BioNetNode *bn_node_tmp = module->GetNodeAttr(i)->clone();
		bn_node_tmp->original_node_id = gcg->GetNodeAttr(sub_node_id)->original_node_id;
		module_to_extension[i] = editor_tmp->InsertNode(bn_node_tmp);
		// Extend the visual information list
		if (visual_info_list != NULL && module_visual_info_list != NULL) {
			double x = module_visual_info_list->at(i).pos_x;
			double y = module_visual_info_list->at(i).pos_y;
			double w = module_visual_info_list->at(i).width;
			double h = module_visual_info_list->at(i).height;
			if (w_c/w_z < h_c/h_z) { // scale by x-axis
				double scale = w_c/w_z;
				(*visual_info_list)[module_to_extension[i]] = new NodeVisualInfo((x-x_min)*scale + x_c, (y-y_min)*scale + y_c + h_c/2 - h_z*scale/2, w*scale, h*scale, scale);
			}
			else {	// scale by y-axis
				double scale = h_c/h_z;
				(*visual_info_list)[module_to_extension[i]] = new NodeVisualInfo((x-x_min)*scale + x_c + w_c/2 - w_z*scale/2 , (y-y_min)*scale + y_c, w*scale, h*scale, scale);
			}
		}
	}
	// Add edges
	for (int i = 0; i < module->NodeCount(); i++)
		for (int j = 0; j < module->NodeCount(); j++)
			if (i != j && module->GetEdgeAttr(i,j) != NULL)
				editor_tmp->InsertEdge(module_to_extension[i], module_to_extension[j], module->GetEdgeAttr(i,j)->clone());
	// TODO: permute inputs and fix for the multiple output case
	// Add connections to the output
	int adj_to_output_node_id = gcg->GetOutEdge(sub_node_id, 0);
	editor_tmp->InsertEdge(module_to_extension[module->outputs[0]], circuit_to_extension[adj_to_output_node_id], gcg->GetEdgeAttr(sub_node_id, adj_to_output_node_id)->clone());
	// Add connections to the inputs
	int number_of_inputs = gcg->InEdgeCount(sub_node_id);
	if (number_of_inputs == 1) {
		int adj_to_input_node_id = gcg->GetInEdge(sub_node_id, 0);
		editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id], module_to_extension[module->inputs[0]], gcg->GetEdgeAttr(adj_to_input_node_id, sub_node_id)->clone());
	}
	else if (number_of_inputs == 2) {
		int adj_to_input_node_id_0 = gcg->GetInEdge(sub_node_id, 0);
		int adj_to_input_node_id_1 = gcg->GetInEdge(sub_node_id, 1);
		if (module->inputs.size() == 1) {	// two signal inputs connect to the same input node in the module, e.g., hrpR --------> hrpR-hrpS <-------- hrpS
			editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id_0], module_to_extension[module->inputs[0]], new BioNetEdge(Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_0)->component, module->GetNodeAttr(module->inputs[0])->component, partDB)));
			editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id_1], module_to_extension[module->inputs[0]], new BioNetEdge(Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_1)->component, module->GetNodeAttr(module->inputs[0])->component, partDB)));
		}
		else {
			if (Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_0)->component, module->GetNodeAttr(module->inputs[0])->component, partDB) != NONE) {
				editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id_0], module_to_extension[module->inputs[0]], new BioNetEdge(Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_0)->component, module->GetNodeAttr(module->inputs[0])->component, partDB)));
				editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id_1], module_to_extension[module->inputs[1]], new BioNetEdge(Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_1)->component, module->GetNodeAttr(module->inputs[1])->component, partDB)));
			}
			else {
				editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id_0], module_to_extension[module->inputs[1]], new BioNetEdge(Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_0)->component, module->GetNodeAttr(module->inputs[1])->component, partDB)));
				editor_tmp->InsertEdge(circuit_to_extension[adj_to_input_node_id_1], module_to_extension[module->inputs[0]], new BioNetEdge(Find_Regulation(gcg->GetNodeAttr(adj_to_input_node_id_1)->component, module->GetNodeAttr(module->inputs[0])->component, partDB)));
			}
		}
	}
	else
		cout << "WE HAVE NOT PROCESSED YET FOR MODULES THAT HAVE MORE THAN TWO INPUTS " << endl;
	return new GeneCircuitGraph(editor_tmp);
}

void Framework::Scalability(int number_of_gates, int number_of_inputs, int number_of_DB_replicates, int module_library_size, int number_of_solutions) {
	bool finish;
	//clock_t begin, end;
	do {
		finish = false;
		//begin = clock();
		// Randomize topology
		dp.gcg = *(GenerateGateNetwork(number_of_gates, number_of_inputs, true));
		dp.gcg.TestPrint(&std::cout);
		break;
		ResetRunningTime();
		//finish = TopologyMatch(number_of_solutions, module_library_size, number_of_DB_replicates);
		//end = clock();
	}
	while (!finish);
	cout << getRunningTime() << "\t";
}

vector<std::pair<GeneCircuitGraph*, SimulationData> >* Framework::Optimize(int mode) {
	vector<GeneCircuitGraph*>* solution_list = ModuleMatch(1,0);
	vector<std::pair<GeneCircuitGraph*, SimulationData> >* tmp = new vector<std::pair<GeneCircuitGraph*, SimulationData> >;
	if (solution_list->empty()) {
		cout << "No solution found!"<< endl;
		return tmp;
	}
	else {
		bool is_zero_input_matrix = true;
		for (int row = 0; row < dp.desired_behavior.input_values.size(); row++) {
			for (int col = 0; col < dp.desired_behavior.input_values[row].size(); col++)
				if (dp.desired_behavior.input_values[row][col] != 0) {
					is_zero_input_matrix = false;
					break;
				}
			if (!is_zero_input_matrix)
				break;
		}
		if (is_zero_input_matrix) { // Module matching only
			dp.desired_behavior.input_values.clear();
			dp.desired_behavior.output_values.clear();
			tmp->push_back(std::make_pair(solution_list->at(0), dp.desired_behavior));
			return tmp;
		}
		else {	// Mutant search here
			cout << "OPTIMIZING ..." << endl;
			switch (mode) {
				case SI: {	// Simulation only
					Module* module_tmp = new Module(solution_list->at(0)->ConvertToLoader(), dp.desired_behavior.inputs, dp.desired_behavior.outputs);
					DynamicGeneCircuit dgc(module_tmp);
					dgc.SimulateSteadyState(partDB, dp.desired_behavior.input_values, dp.desired_behavior.output_values);
					tmp->push_back(std::make_pair(solution_list->at(0), dp.desired_behavior));
					return tmp;
					break;
				}
				case ES: {	// Exhaustive search
					return tmp;
					break;
				}
				case MS: {	// Modular search
					return tmp;
					break;
				}
				case GA: {	// Genetic algorithm search
					return tmp;
					break;
				}
				default:
					return tmp;
			}
		}
	}
}

bool all_visitor(int n, node_id* small_graph_node_list, node_id* large_graph_node_list, void* all_graph_isomorphism) {
	GraphIsomorphism graph_isomorphism;
	graph_isomorphism.reserve(n);
	for (int i = 0; i < n; i++){
		NodePair node_pair_tmp(small_graph_node_list[i],large_graph_node_list[i]);
		graph_isomorphism.push_back(node_pair_tmp);
	}
	((vector<GraphIsomorphism>*) all_graph_isomorphism)->push_back(graph_isomorphism);
	return false;
}
