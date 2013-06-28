/*
 * Framework.h
 *
 *  Created on: Aug 17, 2012
 *      Author: linh
 */

#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_

#include "VFLIB/ull_sub_state.h"
#include "VFLIB/vf2_sub_state.h"
#include "VFLIB/match.h"
#include "VFLIB/argedit.h"
#include "DesignProblem.h"
#include "PartDatabase.h"
#include "GeneCircuitGraph.h"
#include "DynamicGeneCircuit.h"

enum Optimize_Mode {
	SI,	// Simulation only
	ES,	// Exhaustive search
	MS,	// Modular search
	GA	// Genetic algorithm search
};

bool all_visitor(int n, node_id* small_graph_node_list, node_id* large_graph_node_list, void* all_graph_isomorphism);

class Matching_State {
public:
	Matching_State() {
		// Do nothing
	}
	Matching_State(GeneCircuitGraph* gcg_, int module_id_, GraphIsomorphismCollection* matching_list_, int matching_id_, int signal_id_): gcg(gcg_),
		module_id(module_id_),
		matching_list(matching_list_),
		matching_id(matching_id_),
		signal_id(signal_id_) {
		// Do nothing
	}
	~Matching_State() {
		free_var(gcg);
		free_var(matching_list);
	}
	GeneCircuitGraph* gcg;
	int module_id;
	GraphIsomorphismCollection *matching_list;
	int matching_id;
	int signal_id; 			// to manage the cross-talk
};

class BranchingNodeRecord {	// contains information for each branching node
public:
	bool is_contained;
	string signal_name;
};

class OutputNodeRecord { 	// contains information for each node
public:
	void Print();

	int  cost;
	string output;
	StrTree signal_list;
	vector<BranchingNodeRecord> branch_node_info;
	// For tracing back
	int module_id;			// id of the module that can match and create this output
	int matching_id;
	IdList pre_node_record_list;
};
typedef vector<OutputNodeRecord> OutputNodeRecordList;

class Framework {
public:
	Framework();
	Framework(PartDatabase* partDB_);
	Framework(PartDatabase* partDB_, DesignProblem dp_);
	virtual ~Framework();
	vector<GeneCircuitGraph*>* ModuleMatch(int number_of_solution, int synthetic_module_library_size, int number_of_DB_replicates = 1);
	vector<std::pair<GeneCircuitGraph*, SimulationData> >* Optimize(int mode);
	void Scalability(int number_of_gates, int number_of_inputs, int number_of_DB_replicates, int module_library_size, int number_of_solutions);

	double getRunningTime() {
		return running_time;
	}
	void ResetRunningTime() {
		running_time = 0;
	}
	void addRunningTime(double extra_time) {
		running_time += extra_time;
	}

	PartDatabase* partDB;
	DesignProblem dp;

private:
	double running_time;
};

vector<Molecule*>* SrcComponent2Molecule(GeneCircuitComponent *component, PartDatabase* partDB);
Molecule* DestComponent2Molecule(GeneCircuitComponent *component, PartDatabase* partDB);
int Find_Regulation(GeneCircuitComponent* src, GeneCircuitComponent* dest, PartDatabase *partDB);
bool Check_Update(GeneCircuitGraph *gcg, int src_node, int dest_node, PartDatabase *partDB);
GeneCircuitGraph* Layout(GeneCircuitGraph *gcg, PartDatabase* partDB, map<int,NodeVisualInfo*>* visual_info_list = NULL);
GeneCircuitGraph* GeneCircuitGraphExtension(GeneCircuitGraph* gcg, int sub_node_id, Module* module, PartDatabase *partDB, map<int,NodeVisualInfo*>* visual_info_list, vector<NodeVisualInfo> *module_visual_info_list);
void GCG2StrTree(GeneCircuitGraph* gcg, StrTree* str_tree);
bool check_and_merge(StrTree* original_str_tree, StrTree* sub_str_tree, CheckList* skip_list);	// check if there is a cross-talk effect between them and then merge the second tree into the first one
void UnionStrList(StringList* first_str_list, StringList* second_str_list);
void PrintStringTree(StrTree* str_tree);

#endif /* FRAMEWORK_H_ */
