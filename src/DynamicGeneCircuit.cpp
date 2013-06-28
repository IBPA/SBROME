/*
 * DynamicGeneCircuit.cpp
 *
 *  Created on: May 2, 2012
 *      Author: linh
 */

#include <math.h>
#include "DynamicGeneCircuit.h"
#include "GA/GASimpleGA.h"			// we're going to use the simple GA
#include "GA/GA1DArrayGenome.h" 	// and the 2D binary string genome
#include "GA/GABin2DecGenome.h"
#include "METIS/include/metis.h"

#define CHEAP_CUT		1
#define EXPENSIVE_CUT	10

Bound::Bound (double l, double u) {
	lower = l;
	upper = u;
}

Component::Component() {
	sub_component_list.clear();
	inputs.clear();
	outputs.clear();
	linkages.clear();
	variant_id = 1;
}

Component::~Component() {
	// Do nothing
}

void Component::ClearAll() {
	// TODO:
}

void Component::operator << (ostream& out) {
	for (unsigned int i = 0; i < sub_component_list.size(); i++) {
		sub_component_list[i]->operator << (out);
		out << endl;
	}
}

void Component::Linearize(vector<MoleculePool*>& pool_list) {
	if (sub_component_list.empty()) {
		pool_list.push_back((MoleculePool*)this);
	}
	else
		for (unsigned int i = 0; i < sub_component_list.size(); i++)
			sub_component_list[i]->Linearize(pool_list);
}

MoleculePool::MoleculePool() {
	min_concentration = 0;
	max_concentration = 10000;
	concentration = 0;
	min_concentration_list.clear();
	max_concentration_list.clear();
	current_interval_index.clear();
	sub_component_list.clear();
	children.clear();
	parents.clear();
}

MoleculePool::~MoleculePool() {
	// Do nothing
}

void MoleculePool::operator << (ostream& out) {
	molecule->operator << (out);
	out << "\t variant id = \t" << variant_id << "\t sub-circuit-id = " << sub_circuit_id;
}

double MoleculePool::SteadyStateChange(const PartDatabase* pdb) {
	int type = molecule->getType();
	double old_concentration = concentration;
	switch (type) {
		case LIGAND: {
			if (this->parents.size() > 0) { // Ligand is produced by a metabolism process, we assume that the concentration = 1/2 concentration of enzyme
				// TODO: simulate protein -> ligand
				concentration = this->parents[0]->concentration*0.5;
			}
			break;
		}
		case m_RNA: {
			// mRNA is produced by the transcription process at different genes
			concentration = pdb->getPromoterStrength(molecule->name, this->variant_id);
			for (unsigned int i = 0; i < parents.size(); i++) {
				TranscriptionRegulation tr = pdb->getPromoterRegulation(parents[i]->molecule->name,molecule->name);
				if (tr.regulator_type == ACTIVATORY)
					concentration *= (1 - 1/(1 + tr.binding_affinity*pow(parents[i]->concentration, tr.cooperativity)));
				else // PROHIBITORY
					concentration /= (1 + tr.binding_affinity*pow(parents[i]->concentration, tr.cooperativity));
			}
			concentration += pdb->getPromoterBasal(molecule->name, this->variant_id);
			break;
		}
		case PROTEIN: {
			// Unbound protein is produced by the translation process
			concentration = 0;
			bool is_output_protein = false;
			for (unsigned int i = 0; i < parents.size(); i++) {
				if (parents[i]->molecule->getType() == m_RNA || parents[i]->molecule->getType() == RNA_COMPLEX)
					concentration += parents[i]->concentration;
				else if (parents[i]->molecule->getType() == POOL) { // A special case for the output when a protein (e.g. gfp) is produced from a pool
					concentration = parents[i]->concentration;
					is_output_protein = true;
					break;
				}
				else if (parents[i]->molecule->getType() != LIGAND) {
					cout << "error in simulation: PROTEIN contains: " << parents[i]->molecule->getType() << " which is produced neither m_RNA nor mRNA_COMPLEX nor LIGAND" << endl;
					return UNKNOWN;
				}
			}
			if (is_output_protein)
				break;
			concentration /= pdb->getDegradationRate(molecule->name);
			for (unsigned int i = 0; i < parents.size(); i++)
				if (parents[i]->molecule->getType() == LIGAND) {
					LigandRegulation lr = pdb->getLigandRegulation(parents[i]->molecule->name, molecule->name);
					// ligands always repress the protein production
					concentration *= (pow(lr.dissociation, lr.cooperativity)/(pow(lr.dissociation, lr.cooperativity) + pow(parents[i]->concentration, lr.cooperativity)));
				}
			break;
		}
		case RNA_COMPLEX: {
			concentration = 1;
			for (unsigned int i = 0; i < parents.size(); i++)
				concentration *= parents[i]->concentration;
			break;
		}
		case PROTEIN_COMPLEX: {
			concentration = 1;
			for (unsigned int i = 0; i < parents.size(); i++)
				concentration *= parents[i]->concentration;
			break;
		}
		case LIGAND_PROTEIN_COMPLEX: {
			// Bound protein is produced by translation process + binding of inducers
			// TODO: Currently, only one protein + one ligand are considered
			string protein_name, ligand_name;
			int protein_position, ligand_position;
			for (unsigned int i = 0; i < parents.size(); i++) {
				if (parents[i]->molecule->getType() == PROTEIN) {
					protein_name = parents[i]->molecule->name;
					protein_position = i;
				}
				else if (parents[i]->molecule->getType() == LIGAND) {
					ligand_name = parents[i]->molecule->name;
					ligand_position = i;
				}
				else {
					cout << "error in simulation" << endl;
					return UNKNOWN;
				}
			}
			// ligands always activate the ligand-protein complex production
			LigandRegulation lr = pdb->getLigandRegulation(parents[ligand_position]->molecule->name, parents[protein_position]->molecule->name);
			concentration = parents[protein_position]->concentration*(pow(parents[ligand_position]->concentration, lr.cooperativity)/(pow(lr.dissociation, lr.cooperativity) + pow(parents[ligand_position]->concentration, lr.cooperativity)));
			break;
		}
		case POOL: {
			concentration = 0;
			for (unsigned int i = 0; i < parents.size(); i++)
				concentration += parents[i]->concentration;
			break;
		}
		case INPUT: {
			// Do nothing
			break;
		}
		case OUTPUT: {
			concentration = 0;
			for (unsigned int i = 0; i < parents.size(); i++)
				concentration += parents[i]->concentration;
			break;
		}
	}

	return (concentration - old_concentration)*(concentration - old_concentration);
}

double MoleculePool::EstimateSteadyStateBound(const PartDatabase* pdb) {
	int type = molecule->getType();
	double old_min = this->min_concentration;
	double old_max = this->max_concentration;
	switch (type) {
		case LIGAND: {
			if (this->parents.size() > 0) { // Ligand is produced by a metabolism process, we assume that the concentration = 1/2 concentration of enzyme
				// TODO: simulate protein -> ligand
				min_concentration = this->parents[0]->min_concentration*0.5;
				max_concentration = this->parents[0]->max_concentration*0.5;
			}
			break;
		}
		case m_RNA: {
			// mRNA is produced by the transcription process at different genes
			// lower bound
			double current_min, current_max;
			int number_of_mutants = pdb->getNumberOfMutants(this->molecule->name);
			for (int mutant = 1; mutant <= number_of_mutants; mutant++) {
				min_concentration = pdb->getPromoterStrength(molecule->name, mutant);
				for (unsigned int i = 0; i < parents.size(); i++) {
					TranscriptionRegulation tr = pdb->getPromoterRegulation(parents[i]->molecule->name,molecule->name);
					if (tr.regulator_type == ACTIVATORY)
						min_concentration *= (1 - 1/(1 + tr.binding_affinity*pow(parents[i]->min_concentration, tr.cooperativity)));
					else // PROHIBITORY
						min_concentration /= (1 + tr.binding_affinity*pow(parents[i]->max_concentration, tr.cooperativity));
				}
				min_concentration += pdb->getPromoterBasal(molecule->name, mutant);
				// upper bound
				max_concentration = pdb->getPromoterStrength(molecule->name, mutant);
				for (unsigned int i = 0; i < parents.size(); i++) {
					TranscriptionRegulation tr = pdb->getPromoterRegulation(parents[i]->molecule->name,molecule->name);
					if (tr.regulator_type == ACTIVATORY)
						max_concentration *= (1 - 1/(1 + tr.binding_affinity*pow(parents[i]->max_concentration, tr.cooperativity)));
					else // PROHIBITORY
						max_concentration /= (1 + tr.binding_affinity*pow(parents[i]->min_concentration, tr.cooperativity));
				}
				max_concentration += pdb->getPromoterBasal(molecule->name, mutant);
				if (mutant == 1) {
					current_min = min_concentration;
					current_max = max_concentration;
				}
				else {
					if (current_min > min_concentration)
						current_min = min_concentration;
					if (current_max < max_concentration)
						current_max = max_concentration;
				}
			}
			min_concentration = current_min;
			max_concentration = current_max;
			break;
		}
		case PROTEIN: {
			// Unbound protein is produced by the translation process
			min_concentration = max_concentration = 0;
			bool is_output_protein = false;
			for (unsigned int i = 0; i < parents.size(); i++) {
				if (parents[i]->molecule->getType() == m_RNA || parents[i]->molecule->getType() == RNA_COMPLEX) {
					min_concentration += parents[i]->min_concentration;
					max_concentration += parents[i]->max_concentration;
				}
				else if (parents[i]->molecule->getType() == POOL) { // A special case for the output when a protein (e.g. gfp) is produced from a pool
					min_concentration = parents[i]->min_concentration;
					max_concentration = parents[i]->max_concentration;
					is_output_protein = true;
					break;
				}
				else if (parents[i]->molecule->getType() != LIGAND) {
					cout << "error in simulation: PROTEIN contains: " << parents[i]->molecule->getType() << " which is produced neither m_RNA nor mRNA_COMPLEX nor LIGAND" << endl;
					return UNKNOWN;
				}
			}
			if (is_output_protein)
				break;
			min_concentration /= pdb->getDegradationRate(molecule->name);
			max_concentration /= pdb->getDegradationRate(molecule->name);
			for (unsigned int i = 0; i < parents.size(); i++)
				if (parents[i]->molecule->getType() == LIGAND) {
					LigandRegulation lr = pdb->getLigandRegulation(parents[i]->molecule->name, molecule->name);
					// ligands always repress the protein production
					min_concentration *= (pow(lr.dissociation, lr.cooperativity)/(pow(lr.dissociation, lr.cooperativity) + pow(parents[i]->min_concentration, lr.cooperativity)));
					max_concentration *= (pow(lr.dissociation, lr.cooperativity)/(pow(lr.dissociation, lr.cooperativity) + pow(parents[i]->max_concentration, lr.cooperativity)));
				}
			break;
		}
		case RNA_COMPLEX: {
			min_concentration = max_concentration = 1;
			for (unsigned int i = 0; i < parents.size(); i++) {
				min_concentration *= parents[i]->min_concentration;
				max_concentration *= parents[i]->max_concentration;
			}
			break;
		}
		case PROTEIN_COMPLEX: {
			min_concentration = max_concentration = 1;
			for (unsigned int i = 0; i < parents.size(); i++) {
				min_concentration *= parents[i]->min_concentration;
				max_concentration *= parents[i]->max_concentration;
			}
			break;
		}
		case LIGAND_PROTEIN_COMPLEX: {
			// Bound protein is produced by translation process + binding of inducers
			// TODO: Currently, only one protein + one ligand are considered
			string protein_name, ligand_name;
			int protein_position, ligand_position;
			for (unsigned int i = 0; i < parents.size(); i++) {
				if (parents[i]->molecule->getType() == PROTEIN) {
					protein_name = parents[i]->molecule->name;
					protein_position = i;
				}
				else if (parents[i]->molecule->getType() == LIGAND) {
					ligand_name = parents[i]->molecule->name;
					ligand_position = i;
				}
				else {
					cout << "error in simulation" << endl;
					return UNKNOWN;
				}
			}
			// ligands always activate the ligand-protein complex production
			LigandRegulation lr = pdb->getLigandRegulation(parents[ligand_position]->molecule->name, parents[protein_position]->molecule->name);
			min_concentration = parents[protein_position]->min_concentration*(pow(parents[ligand_position]->min_concentration, lr.cooperativity)/(pow(lr.dissociation, lr.cooperativity) + pow(parents[ligand_position]->min_concentration, lr.cooperativity)));
			max_concentration = parents[protein_position]->max_concentration*(pow(parents[ligand_position]->max_concentration, lr.cooperativity)/(pow(lr.dissociation, lr.cooperativity) + pow(parents[ligand_position]->max_concentration, lr.cooperativity)));
			break;
		}
		case POOL: {
			min_concentration = max_concentration = 0;
			for (unsigned int i = 0; i < parents.size(); i++) {
				min_concentration += parents[i]->min_concentration;
				max_concentration += parents[i]->max_concentration;
			}
			break;
		}
	}

	return (min_concentration - old_min)*(min_concentration - old_min) + (max_concentration - old_max)*(max_concentration - old_max);
}

DynamicGeneCircuit::DynamicGeneCircuit() {
	// TODO Auto-generated constructor stub
}

DynamicGeneCircuit::DynamicGeneCircuit(Module* gene_circuit_graph) {
	root_component = new Component;

	int n = gene_circuit_graph->NodeCount();
	root_component->sub_component_list.resize(n);
	for (int i = 0; i < n; i++)
		root_component->sub_component_list[i] = BioNetNode2MoleculePool(gene_circuit_graph->GetNodeAttr(i));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if ((i != j) && gene_circuit_graph->GetEdgeAttr(i,j) != NULL) {
				((MoleculePool*)root_component->sub_component_list[i])->children.push_back((MoleculePool*)root_component->sub_component_list[j]);
				((MoleculePool*)root_component->sub_component_list[j])->parents.push_back((MoleculePool*)root_component->sub_component_list[i]);
			}

	// No linkages since this will return an elementary component
	root_component->linkages.empty();
	// inputs
	for (unsigned int i = 0; i < gene_circuit_graph->inputs.size(); i++)
		root_component->inputs.push_back((MoleculePool*)root_component->sub_component_list[gene_circuit_graph->inputs[i]]);
	// outputs
	for (unsigned int i = 0; i < gene_circuit_graph->outputs.size(); i++)
		root_component->outputs.push_back((MoleculePool*)root_component->sub_component_list[gene_circuit_graph->outputs[i]]);
}

DynamicGeneCircuit::DynamicGeneCircuit(Module* gene_circuit_graph, int number_of_sub_circuits) {
	root_component = new Component;

	int n = gene_circuit_graph->NodeCount();
	root_component->sub_component_list.resize(n);
	for (int i = 0; i < n; i++)
		root_component->sub_component_list[i] = BioNetNode2MoleculePool(gene_circuit_graph->GetNodeAttr(i));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if ((i != j) && gene_circuit_graph->GetEdgeAttr(i,j) != NULL) {
				((MoleculePool*)root_component->sub_component_list[i])->children.push_back((MoleculePool*)root_component->sub_component_list[j]);
				((MoleculePool*)root_component->sub_component_list[j])->parents.push_back((MoleculePool*)root_component->sub_component_list[i]);
			}

	// inputs
	for (unsigned int i = 0; i < gene_circuit_graph->inputs.size(); i++)
		root_component->inputs.push_back((MoleculePool*)root_component->sub_component_list[gene_circuit_graph->inputs[i]]);
	// outputs
	for (unsigned int i = 0; i < gene_circuit_graph->outputs.size(); i++)
		root_component->outputs.push_back((MoleculePool*)root_component->sub_component_list[gene_circuit_graph->outputs[i]]);

	int **adj_matrix = new int*[n];
	for (int i = 0; i < n; i++) {
		adj_matrix[i] = new int[n];
		for (int j = 0; j < n; j++)
			adj_matrix[i][j] = 0;
	}
	int* node_map =  new int[n];
	//int* node_weight = new int[n];
	int number_of_edges = 0;
	for (int i = 0; i < n; i++) {
		node_map[i] = i;
		for (int j = 0; j < n; j++)
			if (i != j && gene_circuit_graph->GetEdgeAttr(i,j) != NULL) {
				// TODO: Should assign weights for edge when we convert from a directed graph -> an undirected one
				//if (gene_circuit_graph->GetEdgeAttr(i,j)->type != NONE)
				if (gene_circuit_graph->GetNodeAttr(i)->original_node_id != gene_circuit_graph->GetNodeAttr(j)->original_node_id) {
					adj_matrix[i][j] = CHEAP_CUT;
					adj_matrix[j][i] = CHEAP_CUT;
				}
				else {
					adj_matrix[i][j] = EXPENSIVE_CUT;
					adj_matrix[j][i] = EXPENSIVE_CUT;
				}
				number_of_edges++;
			}
	}
	// Organize the graph with the input format of METIS
	int* adj_position = new int[n + 1];
	int* adj_array = new int[2*number_of_edges];
	int* adj_weight = new int[2*number_of_edges];
	int current_position = 0;
	for (int i = 0; i < n; i++){
		adj_position[i] = current_position;
		for (int j = 0; j < n; j++)
			if (adj_matrix[i][j] != 0 || adj_matrix[j][i] != 0) {
				//cout << i << " " << j << endl;
				adj_array[current_position] = j;
				adj_weight[current_position] = (adj_matrix[i][j] + adj_matrix[j][i]);
				current_position++;
			}
	}
	adj_position[n] = current_position;

	//for (int i = 0; i < n; i++)
	//	cout << adj_position[i] << " ";
	//cout << endl;
	//for (int i = 0; i < 2*n; i++)
	//	cout << adj_array[i] << " ";
	//cout << endl;
	//for (int i = 0; i < 2*n; i++)
	//	cout << adj_weight[i] << " ";
	//cout << endl;

	int number_of_partitions = number_of_sub_circuits, ncon = 1, obj_value;
	int* partition_id = new int[n];
	for (int i = 0; i < n; i++)
		partition_id[i] = UNKNOWN;

	METIS_PartGraphKway(&n, &ncon, adj_position, adj_array, NULL, NULL, adj_weight, &number_of_partitions, NULL,NULL, NULL, &obj_value, partition_id);
	//cout << endl << obj_value << endl;
	// Justify the node partition
	for (int i = 0; i < n; i++) {
		IdList adj_tmp;
		for (int j = 0; j < n; j++)
			if (i != j && adj_matrix[i][j] == EXPENSIVE_CUT && partition_id[i] != partition_id[j])
				adj_tmp.push_back(j);
		bool is_updated = false;
		for (unsigned int k = 0; k < adj_tmp.size(); k++) {
			for (unsigned int l = k + 1; l < adj_tmp.size(); l++)
				if (partition_id[adj_tmp[k]] == partition_id[adj_tmp[l]]) {
					partition_id[i] = partition_id[adj_tmp[k]];
					is_updated = true;
					break;
				}
			if (is_updated)
				break;
		}
	}
	for (int i = 0; i < n; i++)
		root_component->sub_component_list[i]->sub_circuit_id = partition_id[i];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			if (i != j && gene_circuit_graph->GetEdgeAttr(i,j) != NULL && partition_id[i] != partition_id[j]) {
				bool is_input = false;
				for (unsigned int l = 0; l < gene_circuit_graph->inputs.size(); l++)
					if (i == gene_circuit_graph->inputs[l]) {
						is_input = true;
						break;
					}
				if (!is_input)
					root_component->linkages.push_back((MoleculePool*)root_component->sub_component_list[i]);
			}
	}
	this->PrintOut(&std::cout);
	delete []partition_id;
	// TODO: garbage
}

DynamicGeneCircuit::~DynamicGeneCircuit() {
	RecursiveClear(root_component);
}

void DynamicGeneCircuit::RecursiveClear(Component* component_node) {
	for(unsigned int i = 0; i < component_node->sub_component_list.size(); i++)
		RecursiveClear(component_node->sub_component_list[i]);
	delete component_node;
}

void DynamicGeneCircuit::PrintOut(ostream* f) {
	root_component->operator <<(*f);
}

SimulationStatus DynamicGeneCircuit::SimulateSteadyState(const PartDatabase* pdb, const ValueMatrix& input_values, ValueMatrix& output_values, double threshold_value, int max_iteration) {

	// Resize the output value matrix
	output_values.resize(input_values.size());
	for (unsigned int i = 0; i < output_values.size(); i++)
		output_values[i].resize(root_component->outputs.size());

	bool over_time = false;
	vector<MoleculePool*> molecule_list;
	root_component->Linearize(molecule_list);
	int n = molecule_list.size();
	for (int i = 0; i < n; i++)
		molecule_list[i]->best_concentration_list.clear();

	for (unsigned int condition = 0; condition < input_values.size(); condition++) {
		// Assign the input concentration for inducers
		for (unsigned int i = 0; i < root_component->inputs.size(); i++)
			root_component->inputs[i]->concentration = input_values[condition][i];

		// Simulate here
		int count = 0;
		double error_tmp_max;
		do {
			error_tmp_max = 0;
			for (int i = 0; i < n; i++) {
				double diff = molecule_list[i]->SteadyStateChange(pdb);
				if (diff > error_tmp_max)
					error_tmp_max = diff;
			}
			count++;
		}
		while (error_tmp_max > threshold_value && count < max_iteration);
		if (count >= max_iteration) {
			cout << "Simulation runs over time" << endl;
			over_time = true;
		}
		// debug
		for (int i = 0; i < n; i++)
			cout << molecule_list[i]->molecule->name << "\t\t" << molecule_list[i]->concentration << endl;
		cout << "---------------" << endl;
		// end debug
		for (unsigned int i = 0; i < root_component->outputs.size(); i++)
			output_values[condition][i] = root_component->outputs[i]->concentration;
	}

	return (over_time)?OVER_TIME:CONVERGENCE;
}

SimulationStatus DynamicGeneCircuit::PartialSimulateSteadyState(const PartDatabase* pdb, int sub_circuit_id, double threshold_value, int max_iteration) {
	bool over_time = false;
	vector<MoleculePool*> molecule_list;
	root_component->Linearize(molecule_list);
	int n = molecule_list.size();

	// Simulate here
	int count = 0;
	double error_tmp_max;
	do {
		error_tmp_max = 0;
		for (int i = 0; i < n; i++)
			if (molecule_list[i]->sub_circuit_id == sub_circuit_id) {
				double diff = molecule_list[i]->SteadyStateChange(pdb);
				if (diff > error_tmp_max)
					error_tmp_max = diff;
			}
		count++;
	}
	while (error_tmp_max > threshold_value && count < max_iteration);
	if (count >= max_iteration) {
		cout << "Simulation runs over time" << endl;
		over_time = true;
	}
	return (over_time)?OVER_TIME:CONVERGENCE;
}
SimulationStatus DynamicGeneCircuit::EstimateSteadyStateBound(const PartDatabase* pdb, const ValueMatrix& input_values, ValueMatrix& linkage_lb, ValueMatrix& linkage_ub, double threshold_value, int max_iteration) {

	// Resize the output value matrix
	linkage_lb.resize(input_values.size());
	linkage_ub.resize(input_values.size());
	for (unsigned int i = 0; i < linkage_lb.size(); i++) {
		linkage_lb[i].resize(root_component->linkages.size());
		linkage_ub[i].resize(root_component->linkages.size());
	}

	bool over_time = false;
	vector<MoleculePool*> molecule_list;
	root_component->Linearize(molecule_list);
	int n = molecule_list.size();
	for (int i = 0; i < n; i++)
		molecule_list[i]->best_concentration_list.clear();

	for (unsigned int condition = 0; condition < input_values.size(); condition++) {
		// Assign the input concentration for inducers
		for (unsigned int i = 0; i < root_component->inputs.size(); i++)
			root_component->inputs[i]->min_concentration = root_component->inputs[i]->max_concentration = input_values[condition][i];

		// Simulate here
		int count = 0;
		double error_tmp_max;
		do {
			error_tmp_max = 0;
			for (int i = 0; i < n; i++) {
				double diff = molecule_list[i]->EstimateSteadyStateBound(pdb);
				if (diff > error_tmp_max)
					error_tmp_max = diff;
			}
			count++;
		}
		while (error_tmp_max > threshold_value && count < max_iteration);
		if (count >= max_iteration) {
			cout << "Simulation runs over time" << endl;
			over_time = true;
		}
		// debug
		//for (int i = 0; i < n; i++)
		//	cout << molecule_list[i]->molecule->name << "\t\t" << molecule_list[i]->concentration << endl;
		//cout << "---------------" << endl;
		// end debug
		for (unsigned int i = 0; i < root_component->linkages.size(); i++) {
			linkage_lb[condition][i] = root_component->linkages[i]->min_concentration;
			linkage_ub[condition][i] = root_component->linkages[i]->max_concentration;
		}
		//for (int i = 0; i < root_component->sub_component_list.size(); i++) {
		//	cout << ((MoleculePool*)root_component->sub_component_list[i])->molecule->name << "\t" << ((MoleculePool*)root_component->sub_component_list[i])->min_concentration << "\t" << ((MoleculePool*)root_component->sub_component_list[i])->max_concentration << endl;
		//}
	}

	return (over_time)?OVER_TIME:CONVERGENCE;
}

DynamicGeneCircuit* DummyForGA::dynamic_gene_circuit = NULL;
const PartDatabase* DummyForGA::part_database = NULL;
//const DesignProblem* DummyForGA::design_problem = NULL;
const ValueMatrix* DummyForGA::inputs = new ValueMatrix;
const ValueMatrix* DummyForGA::outputs = new ValueMatrix;
vector<MoleculePool*> DummyForGA::promoter_list;
IdList DummyForGA::number_of_mutants_list;

void DynamicGeneCircuit::ExhaustiveSearch(const PartDatabase* pdb, const ValueMatrix& input_values, const ValueMatrix& output_values) {
	// Find all promoters that have more than one mutant

	vector<MoleculePool*> molecule_list;
	root_component->Linearize(molecule_list);
	DummyForGA::promoter_list.clear();
	DummyForGA::number_of_mutants_list.clear();
	for (unsigned int i = 0; i < molecule_list.size(); i++)
		if (molecule_list[i]->molecule->getType() == m_RNA) {
			int number_of_mutants = pdb->getNumberOfMutants(molecule_list[i]->molecule->name);
			if (number_of_mutants > 1) {
				DummyForGA::promoter_list.push_back(molecule_list[i]);
				DummyForGA::number_of_mutants_list.push_back(number_of_mutants);
			}
			// TEST HERE
			//cout << molecule_list[i]->molecule->name << "\t" << number_of_mutants << endl;
		}
	DummyForGA::dynamic_gene_circuit = this;
	DummyForGA::part_database = pdb;
	DummyForGA::inputs = &input_values;

	// Brute-force search
	IdList best_mutant (DummyForGA::promoter_list.size());
	for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
		(DummyForGA::promoter_list[i])->variant_id = 1;

	bool terminate = false;
	double error = 1e10;
	int count = 0;
	do {
		// TEST HERE
		count++;
		if (count % 1000 == 0)
			cout << count << endl;

		double error_tmp = 0;
		bool good_simulation = true;

		ValueMatrix output_value_tmp;
		good_simulation = good_simulation && (DummyForGA::dynamic_gene_circuit->SimulateSteadyState(DummyForGA::part_database, input_values, output_value_tmp) == CONVERGENCE);
		for (unsigned int i = 0; i < output_value_tmp.size(); i++)
			for (unsigned int j = 0; j < output_value_tmp[i].size(); j++)
				error_tmp += pow(output_value_tmp[i][j] - output_values[i][j],2);
		// Estimate error
		if (good_simulation && error_tmp < error) {
			// Update the best solution
			//cout << error_tmp << endl;
			error = error_tmp;
			for (unsigned int i = 0; i < best_mutant.size(); i++)
				best_mutant[i] = (DummyForGA::promoter_list[i])->variant_id;
		}
		terminate = true;
		for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
			if ((DummyForGA::promoter_list[i])->variant_id < DummyForGA::number_of_mutants_list[i]) {
				(DummyForGA::promoter_list[i])->variant_id++;
				terminate = false;
				break;
			}
			else
				(DummyForGA::promoter_list[i])->variant_id = 1;
	}
	while (!terminate);

	// Simulate again
	for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
		(DummyForGA::promoter_list[i])->variant_id = best_mutant[i];
	ValueMatrix output_values_tmp;
	SimulateSteadyState(DummyForGA::part_database, *(DummyForGA::inputs), output_values_tmp);
	for (unsigned int i = 0; i < output_values_tmp.size(); i++) {
		for (unsigned int j = 0; j < output_values_tmp[i].size(); j++)
			cout << output_values_tmp[i][j] << "(" << output_values[i][j] << ")" << "\t";
		cout << endl;
	}
}

void DynamicGeneCircuit::ModularSearch(const PartDatabase* pdb, const ValueMatrix& input_values, const ValueMatrix& output_values) {
	// Estimate bounds
	ValueMatrix linkage_lb, linkage_ub;
	EstimateSteadyStateBound(pdb, input_values, linkage_lb, linkage_ub);
	int number_of_sub_circuits = 2;
	int resolution = 5;
	ValueList best_linkage_value_list(input_values.size(), 0);

	for (unsigned int condition = 0; condition < input_values.size(); condition++) {	// For each condition
		cout << "------------------------" << endl;
		// Assume that there is only one linkage and two-subcircuits
		if (root_component->linkages.size() != 1) {
			cout << "THERE ARE " << root_component->linkages.size() << ", NOT CONSIDER YET!.";
			return;
		}
		double lb = linkage_lb[condition][0], ub = linkage_ub[condition][0];
		cout << "lb = " << lb << "---- ub = " << ub << endl;
		double best_global_error = 1e10;
		for (int r = 0; r < resolution; r++) {
			double desired_linkage_value = lb + (ub - lb)*r/((double)(resolution - 1));
			// Update the global inputs
			for (unsigned int i = 0; i < input_values[condition].size(); i++)
				root_component->inputs[i]->concentration = input_values[condition][i];
			//for (int sub_circuit = number_of_sub_circuits - 1; sub_circuit >= 0; sub_circuit--) {
			for (int sub_circuit = 0; sub_circuit < number_of_sub_circuits; sub_circuit++) {
				vector<MoleculePool*> molecule_list;
				root_component->Linearize(molecule_list);
				DummyForGA::promoter_list.clear();
				DummyForGA::number_of_mutants_list.clear();
				for (unsigned int i = 0; i < molecule_list.size(); i++)
					if (molecule_list[i]->molecule->getType() == m_RNA && molecule_list[i]->sub_circuit_id == sub_circuit) {
						int number_of_mutants = pdb->getNumberOfMutants(molecule_list[i]->molecule->name);
						if (number_of_mutants > 1) {
							DummyForGA::promoter_list.push_back(molecule_list[i]);
							DummyForGA::number_of_mutants_list.push_back(number_of_mutants);
						}
						// TEST HERE
						//cout << molecule_list[i]->molecule->name << "\t" << number_of_mutants << endl;
					}
				DummyForGA::dynamic_gene_circuit = this;
				DummyForGA::part_database = pdb;
				DummyForGA::inputs = &input_values;

				for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
					(DummyForGA::promoter_list[i])->variant_id = 1;

				int count = 0;
				bool terminate = false;
				double error = 1e10;
				double output_tmp;
				double best_linkage_value;
				do {
					// TEST HERE
					count++;
					if (count % 1000 == 0)
						cout << count << endl;

					double error_tmp = 0;
					bool good_simulation = true;

					good_simulation = good_simulation && (DummyForGA::dynamic_gene_circuit->PartialSimulateSteadyState(DummyForGA::part_database, sub_circuit) == CONVERGENCE);

					//if (condition == 2 && sub_circuit == 0) {
					//	cout << "------------------------------\n----------------------------------" << endl;
					//	for (int i = 0; i < root_component->sub_component_list.size(); i++)
					//		cout << ((MoleculePool*)root_component->sub_component_list[i])->molecule->name << "\t" << ((MoleculePool*)root_component->sub_component_list[i])->concentration << endl;
					//}

					//if (sub_circuit == 1)
					if (sub_circuit == 0)
						error_tmp += pow(root_component->linkages[0]->concentration - desired_linkage_value,2);
					else
						error_tmp += pow(root_component->outputs[0]->concentration - output_values[condition][0],2);
					// Estimate error
					if (good_simulation && error_tmp < error) {
						// Update the best solution
						//cout << error_tmp << endl;
						error = error_tmp;
						output_tmp = root_component->outputs[0]->concentration;
						if (sub_circuit == 0) //if (sub_circuit == 1)
							best_linkage_value = root_component->linkages[0]->concentration;
					}
					terminate = true;
					for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
						if ((DummyForGA::promoter_list[i])->variant_id < DummyForGA::number_of_mutants_list[i]) {
							(DummyForGA::promoter_list[i])->variant_id++;
							terminate = false;
							break;
						}
						else
							(DummyForGA::promoter_list[i])->variant_id = 1;
				}
				while (!terminate);
				// Update linkage value
				if (sub_circuit == 0) //if (sub_circuit == 1)
					root_component->linkages[0]->concentration = best_linkage_value;
				else
					if (error < best_global_error) {
						best_global_error = error;
						cout << output_tmp << endl;
						// Update the best linkage value
						best_linkage_value_list[condition] = root_component->linkages[0]->concentration;
					}
			}
		}
	}
	// Search best mutants for the combination
	IdMatrix best_mutants;
	best_mutants.resize(number_of_sub_circuits);
	//for (int sub_circuit = number_of_sub_circuits - 1; sub_circuit >= 0; sub_circuit--) {
	for (int sub_circuit = 0; sub_circuit < number_of_sub_circuits; sub_circuit++) {
		vector<MoleculePool*> molecule_list;
		root_component->Linearize(molecule_list);
		DummyForGA::promoter_list.clear();
		DummyForGA::number_of_mutants_list.clear();
		for (unsigned int i = 0; i < molecule_list.size(); i++)
			if (molecule_list[i]->molecule->getType() == m_RNA && molecule_list[i]->sub_circuit_id == sub_circuit) {
				int number_of_mutants = pdb->getNumberOfMutants(molecule_list[i]->molecule->name);
				if (number_of_mutants > 1) {
					DummyForGA::promoter_list.push_back(molecule_list[i]);
					DummyForGA::number_of_mutants_list.push_back(number_of_mutants);
				}
				// TEST HERE
				//cout << molecule_list[i]->molecule->name << "\t" << number_of_mutants << endl;
			}
		best_mutants[sub_circuit].resize(DummyForGA::promoter_list.size());
		int count = 0;
		bool terminate = false;
		double error = 1e10;
		//double best_linkage_value;
		ValueList best_linkage_value_list_tmp(input_values.size(), 0);

		for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
			(DummyForGA::promoter_list[i])->variant_id = 1;
		do {
			double error_tmp = 0;
			bool good_simulation = true;
			ValueList linkage_value_list_tmp(input_values.size(), 0);
			for (unsigned int condition = 0; condition < input_values.size(); condition++) {	// For each condition
				double best_global_error = 1e10;
				for (unsigned int i = 0; i < input_values[condition].size(); i++)
					root_component->inputs[i]->concentration = input_values[condition][i];
				root_component->linkages[0]->concentration = best_linkage_value_list_tmp[condition];
				// TEST HERE
				count++;
				if (count % 100 == 0)
					cout << count << endl;

				good_simulation = good_simulation && (DummyForGA::dynamic_gene_circuit->PartialSimulateSteadyState(DummyForGA::part_database, sub_circuit) == CONVERGENCE);
				if (sub_circuit == 0) { //if (sub_circuit == 1) {
					error_tmp += pow(root_component->linkages[0]->concentration - best_linkage_value_list[condition],2);
					linkage_value_list_tmp[condition] = root_component->linkages[0]->concentration;
				}
				else
					error_tmp += pow(root_component->outputs[0]->concentration - output_values[condition][0],2);
			}
			// Estimate error
			if (good_simulation && error_tmp < error) {
				// Update the best solution
				//cout << error_tmp << endl;
				error = error_tmp;
				// Update best mutants
				for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
					best_mutants[sub_circuit][i] = (DummyForGA::promoter_list[i])->variant_id;
				if (sub_circuit == 0) //if (sub_circuit == 1)
					for (unsigned int condition = 0; condition < input_values.size(); condition++)
						best_linkage_value_list_tmp[condition] = linkage_value_list_tmp[condition];
			}
			terminate = true;
			for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
				if ((DummyForGA::promoter_list[i])->variant_id < DummyForGA::number_of_mutants_list[i]) {
					(DummyForGA::promoter_list[i])->variant_id++;
					terminate = false;
					break;
				}
				else
					(DummyForGA::promoter_list[i])->variant_id = 1;
		}
		while (!terminate);
	}

	vector<MoleculePool*> molecule_list;
	root_component->Linearize(molecule_list);
	DummyForGA::promoter_list.clear();
	DummyForGA::number_of_mutants_list.clear();
	for (int sub_circuit = 0; sub_circuit < number_of_sub_circuits; sub_circuit++) {
		int tmp = -1;
		for (unsigned int i = 0; i < molecule_list.size(); i++)
			if (molecule_list[i]->molecule->getType() == m_RNA && molecule_list[i]->sub_circuit_id == sub_circuit) {
				int number_of_mutants = pdb->getNumberOfMutants(molecule_list[i]->molecule->name);
				if (number_of_mutants > 1) {
					DummyForGA::promoter_list.push_back(molecule_list[i]);
					tmp++;
					DummyForGA::promoter_list[DummyForGA::promoter_list.size() - 1]->variant_id = best_mutants[sub_circuit][tmp];
				}
		}
	}
	DummyForGA::dynamic_gene_circuit = this;
	DummyForGA::part_database = pdb;
	DummyForGA::inputs = &input_values;
	ValueMatrix output_values_tmp;
	SimulateSteadyState(DummyForGA::part_database, *(DummyForGA::inputs), output_values_tmp);
	for (unsigned int i = 0; i < output_values_tmp.size(); i++) {
		for (unsigned int j = 0; j < output_values_tmp[i].size(); j++)
			cout << output_values_tmp[i][j]  << "\t";
		cout << endl;
	}
}

void DynamicGeneCircuit::GASearch(const PartDatabase* pdb, const ValueMatrix& input_values, const ValueMatrix& output_values) {
	// Find all promoters that have more than one mutant

	vector<MoleculePool*> molecule_list;
	root_component->Linearize(molecule_list);
	DummyForGA::promoter_list.clear();
	DummyForGA::number_of_mutants_list.clear();
	for (unsigned int i = 0; i < molecule_list.size(); i++)
		if (molecule_list[i]->molecule->getType() == m_RNA) {
			int number_of_mutants = pdb->getNumberOfMutants(molecule_list[i]->molecule->name);
			if (number_of_mutants > 1) {
				DummyForGA::promoter_list.push_back(molecule_list[i]);
				DummyForGA::number_of_mutants_list.push_back(number_of_mutants);
			}
		}
	// Run GA here
	// Number of genes for GA
	//int length = 0;
	GABin2DecPhenotype map;
	for (unsigned int i = 0; i < DummyForGA::number_of_mutants_list.size(); i++)
		map.add(ceil(log2(DummyForGA::number_of_mutants_list[i])), 0, DummyForGA::number_of_mutants_list[i] - 1);

	DummyForGA::dynamic_gene_circuit = this;
	DummyForGA::part_database = pdb;
	//DummyForGA::design_problem = &design_problem;
	DummyForGA::inputs = &input_values;
	DummyForGA::outputs = &output_values;

	//GA1DArrayGenome<unsigned int> genome(length, DummyForGA::GAFitness);				// create a genome
	GABin2DecGenome genome(map, DummyForGA::GAFitness);
	GASimpleGA ga(genome);																// create the genetic algorithm
	ga.populationSize(20); // 200
	ga.nGenerations(5);	// 50
	ga.pMutation(0.01);
	ga.pCrossover(0.9);

	ga.evolve();
	//cout << ga.statistics() << endl;

	// Simulate again
	GAGenome g = ga.statistics().bestIndividual();
	//cout << ga.statistics().bestIndividual();
	for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
		DummyForGA::promoter_list[i]->variant_id = ((GA1DArrayGenome<unsigned int>*)&g)->gene(i) % DummyForGA::number_of_mutants_list[i] + 1;
	ValueMatrix output_values_tmp;
	SimulateSteadyState(DummyForGA::part_database, *(DummyForGA::inputs), output_values_tmp);
	for (unsigned int i = 0; i < output_values_tmp.size(); i++) {
		for (unsigned int j = 0; j < output_values_tmp[i].size(); j++)
			cout << output_values_tmp[i][j]  << "\t";
		cout << endl;
	}
}

float DummyForGA::GAFitness(GAGenome& g) {
	GA1DArrayGenome<unsigned int> &genome = (GA1DArrayGenome<unsigned int> &) g;
	for (unsigned int i = 0; i < DummyForGA::promoter_list.size(); i++)
		(DummyForGA::promoter_list[i])->variant_id = genome.gene(i) % DummyForGA::number_of_mutants_list[i] + 1;

	ValueMatrix output_values;
	DummyForGA::dynamic_gene_circuit->SimulateSteadyState(DummyForGA::part_database, *(DummyForGA::inputs), output_values);
	double error = 0;
	for (unsigned int i = 0; i < output_values.size(); i++)
		for (unsigned int j = 0; j < output_values[i].size(); j++)
			error += pow(output_values[i][j] - (*DummyForGA::outputs)[i][j],2);
	//cout << 100 - error << endl;
	return 100 - error;
}

MoleculePool* BioNetNode2MoleculePool(BioNetNode* bio_node) {
	MoleculePool* tmp = new MoleculePool;
	tmp->molecule = (Molecule*)bio_node->component->clone();
	return tmp;
}
