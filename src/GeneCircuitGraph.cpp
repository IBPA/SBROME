/*
 * GeneCircuitGraph.cpp
 *
 *  Created on: May 2, 2012
 *      Author: linh, UC Davis
 */

#include <stdlib.h>
#include <stdio.h>
#include "GeneCircuitGraph.h"

string Id2Str(int type_id) {
	switch (type_id) {
		case LIGAND:
			return "Ligand";
		case m_RNA:
			return "mRNA";
		case t_RNA:
			return "tRNA";
		case PROTEIN:
			return "Protein";
		case POOL:
			return "Pool";
		case LIGAND_PROTEIN_COMPLEX:
			return "LigandProteinComplex";
		case RNA_COMPLEX:
			return "RNAComplex";
		case PROTEIN_COMPLEX:
			return "ProteinComplex";
		case INPUT:
			return "INPUT";
		case OUTPUT:
			return "OUTPUT";
		case LINKER:
			return "LINKER";
		case YES_GATE:
			return "YES_GATE";
		case NOT_GATE:
			return "NOT_GATE";
		case AND_GATE:
			return "AND_GATE";
		case NAND_GATE:
			return "NAND_GATE";
		case OR_GATE:
			return "OR_GATE";
		case NOR_GATE:
			return "NOR_GATE";
		default:
			cout << "Error: There is no id type: " << type_id << endl;
			return "UNKNOWN";
	}
}

int Str2Id(string type_str) {
	if (type_str.compare("Ligand") == STR_EQ)
		return LIGAND;
	else if (type_str.compare("mRNA") == STR_EQ)
		return m_RNA;
	else if (type_str.compare("mRNA(t)") == STR_EQ)
		return m_RNA;
	else if (type_str.compare("tRNA") == STR_EQ)
			return t_RNA;
	else if (type_str.compare("Protein") == STR_EQ)
		return PROTEIN;
	else if (type_str.compare("Pool") == STR_EQ)
		return POOL;
	else if (type_str.compare("LigandProteinComplex") == STR_EQ)
		return LIGAND_PROTEIN_COMPLEX;
	else {
		cout << "Error: There is no type: " << type_str << endl;
		return UNKNOWN;
	}
}

bool SingleMolecule::CompareTo(GeneCircuitComponent* another_component) {
	if (another_component->getCategory() == MOLECULE) {
		if (another_component->getType() != UNKNOWN_MOLECULE && this->getType() != UNKNOWN_MOLECULE && (another_component->getType() != this->getType()
			|| (another_component->getName().compare("") != STR_EQ && this->getName().compare("") != STR_EQ && another_component->getName().compare(this->getName()) != STR_EQ)))
				return false;
		else
			return true;
	}
	else if (another_component->getCategory() == CIRCUIT_SIGNAL) {
		GeneCircuitComponent *tmp = ((InputSignal*)another_component)->physical_instance;
		if (tmp->getCategory() != MOLECULE
			|| (this->getType() != UNKNOWN_MOLECULE && (tmp->getType() != this->getType()
			|| (tmp->getName().compare("") != STR_EQ && this->getName().compare("") != STR_EQ && tmp->getName().compare(this->getName()) != STR_EQ))))
				return false;
		else
			return true;
	}
	else {	// GATE
		return false;
	}
}

GeneCircuitGraph::GeneCircuitGraph() {
	// Do nothing
}
GeneCircuitGraph::GeneCircuitGraph(ARGLoader* loader): ARGraph<BioNetNode, BioNetEdge>(loader) {
	this->SetNodeDestroyer(new BioNetNodeDestroyer());
	this->SetNodeComparator(new BioNetNodeComparator());
	this->SetEdgeDestroyer(new BioNetEdgeDestroyer());
	this->SetEdgeComparator(new BioNetEdgeComparator());
}
GeneCircuitGraph::GeneCircuitGraph(istream* f){
	Init(f);
}
void GeneCircuitGraph::Init(istream* f) {
	NewAllocator<BioNetNode> node_allocator;
	NewAllocator<BioNetEdge> edge_allocator;
	StreamARGLoader<BioNetNode, BioNetEdge> loader_temp(&node_allocator, &edge_allocator, *f);
	this->setLoader(&loader_temp);
	this->SetNodeDestroyer(new BioNetNodeDestroyer());
	this->SetNodeComparator(new BioNetNodeComparator());
	this->SetEdgeDestroyer(new BioNetEdgeDestroyer());
	this->SetEdgeComparator(new BioNetEdgeComparator());
	number_of_inputs = number_of_outputs = 0;
	for (int i = 0; i < loader_temp.NodeCount(); i++) {
		if (this->InEdgeCount(i) == 0)
			number_of_inputs++;
		if (this->OutEdgeCount(i) == 0)
			number_of_outputs++;
	}
}
GeneCircuitGraph::GeneCircuitGraph(const GeneCircuitGraph& another_gene_circuit_graph): ARGraph<BioNetNode, BioNetEdge>(another_gene_circuit_graph) {
	// Do nothing
}
GeneCircuitGraph& GeneCircuitGraph::operator = (const GeneCircuitGraph& another_gene_circuit_graph) {
	ARGraph<BioNetNode, BioNetEdge>::operator = (another_gene_circuit_graph);
	return *this;
}
GeneCircuitGraph::~GeneCircuitGraph() {
	// Do nothing
}
void GeneCircuitGraph::TestPrint(ostream* f) {
	StreamARGLoader<BioNetNode, BioNetEdge>::write(*f,*this);
}
bool GeneCircuitGraph::IsComplete() {
	for (int i = 0; i < (this->n); i++)
		if (this->GetNodeAttr(i)->component->getCategory() != CIRCUIT_SIGNAL || this->GetNodeAttr(i)->component->getType() != OUTPUT) {
			if (this->GetNodeAttr(i)->component->name.empty())
				return false;
		}
		else {

		}
	return true;
}

int GeneCircuitGraph::getNumberOfInputs() {
	return number_of_inputs;
}
int GeneCircuitGraph::getNumberOfOutputs() {
	return number_of_outputs;
}

int GeneCircuitGraph::getNumberOfModuleNodes() {
	int count_tmp = 0;
	for (int i = 0; i < this->NodeCount(); i++)
		if (this->GetNodeAttr(i)->component->getCategory() == LOGIC_GATE)
			count_tmp++;
	return count_tmp;
}

ARGLoader* GeneCircuitGraph::ConvertToLoader() {
	ARGEdit* tmp = new ARGEdit;
	for (int i = 0; i < this->NodeCount(); i++)
		tmp->InsertNode(this->GetNodeAttr(i)->clone());
	for (int i = 0; i < this->NodeCount(); i++)
		for (int j = 0; j < this->NodeCount(); j++)
			if (this->GetEdgeAttr(i,j) != NULL)
				tmp->InsertEdge(i,j,this->GetEdgeAttr(i,j)->clone());
	return tmp;
}

int GeneCircuitGraph::getNumberOfUnkownNodes() {
	int count = 0;
	for (int i = 0; i < (this->n); i++)
		if (this->GetNodeAttr(i)->component->name.compare("") == STR_EQ)
			count++;
	return count;
}

GeneCircuitGraph* GeneCircuitGraph::clone() {
	return new GeneCircuitGraph(this->ConvertToLoader());
}

BioNetNode* CreateBioNetNode (string type) {
	BioNetNode* tmp = new BioNetNode;
	tmp->component = CreateComponent(type);
	return tmp;
}

BioNetNode* CreateBioNetNode (string type, string name) {
	BioNetNode* tmp = new BioNetNode;
	tmp->component = CreateComponent(type, name);
	return tmp;
}

GeneCircuitComponent* CreateComponent (string type) {
	return CreateComponent(type, "");
}

GeneCircuitComponent* CreateComponent (string type, string name) {
	GeneCircuitComponent *component;
	if (type.compare("Ligand") == STR_EQ || type.compare("mRNA") == STR_EQ || type.compare("tRNA") == STR_EQ || type.compare("Protein") == STR_EQ || type.compare("Pool") == STR_EQ)
		component = new SingleMolecule(Str2Id(type), name);
	else if (type.compare("mRNA(t)") == STR_EQ)
		component = new SingleMolecule(Str2Id(type), name, TANDEM);
	else if (type.compare("LigandProteinComplex") == STR_EQ) {
		int pos = name.find_first_of('-');
		SingleMolecule* m1 = new SingleMolecule(LIGAND, name.substr(0,pos));
		SingleMolecule* m2 = new SingleMolecule(PROTEIN, name.substr(pos + 1));
		component = new ComplexMolecule(m1,m2);
		delete m1;
		delete m2;
	}
	else if (type.compare("ProteinComplex") == STR_EQ) {
		int pos = name.find_first_of('-');
		SingleMolecule* m1 = new SingleMolecule(PROTEIN, name.substr(0,pos));
		SingleMolecule* m2 = new SingleMolecule(PROTEIN, name.substr(pos + 1));
		component = new ComplexMolecule(m1,m2);
		delete m1;
		delete m2;
	}
	else if (type.compare("RNAComplex") == STR_EQ) {
		int pos = name.find_first_of('-');
		SingleMolecule* m1 = new SingleMolecule(m_RNA, name.substr(0,pos));
		SingleMolecule* m2 = new SingleMolecule(m_RNA, name.substr(pos + 1));
		component = new ComplexMolecule(m1,m2);
		delete m1;
		delete m2;
	} else if (type.compare("INPUT") == STR_EQ)
		component = new InputSignal;
	else if (type.compare("OUTPUT") == STR_EQ)
		component = new OutputSignal;
	else if (type.compare("LINKER") == STR_EQ)
			component = new LinkSignal;
	else if (type.compare("YES_GATE") == STR_EQ)
		component = new YESGate;
	else if (type.compare("NOT_GATE") == STR_EQ)
		component = new NOTGate;
	else if (type.compare("AND_GATE") == STR_EQ)
		component = new ANDGate;
	else if (type.compare("NAND_GATE") == STR_EQ)
		component = new NANDGate;
	else if (type.compare("OR_GATE") == STR_EQ)
		component = new ORGate;
	else if (type.compare("NOR_GATE") == STR_EQ)
		component = new NORGate;
	else cout << "Error: " << type << ", no such type defined" << endl;
	component->name = name;
	return component;
}

inline istream& operator >> (istream& in, BioNetNode& node) {
	string type;
	in >> type;
	node.component = CreateComponent(type);
	node.component->operator >>(in);
	return in;
}
inline ostream& operator << (ostream& out, BioNetNode& node) {
	//out << "\t" << node.original_node_id << "\t";
	out << "\t";
	node.component->operator <<(out);
	if (node.component->getType() == m_RNA && node.component->variant_id != UNKNOWN)
		out << "\t" << "variant id = " << node.component->variant_id;
	return out;
}

inline istream& operator >> (istream& in, BioNetEdge& BioNet_edge) {
	string tmp;
	in >> tmp;
	if (tmp.compare("ACT") == 0)
		BioNet_edge.type = ACTIVATORY;
	else if (tmp.compare("REP") == 0)
		BioNet_edge.type = INHIBITORY;
	else if (tmp.compare("UNKNOWN") == 0 || tmp.compare("NA") == 0)
		BioNet_edge.type = UNKNOWN;
	return in;
}

inline ostream& operator << (ostream& out, BioNetEdge BioNet_edge) {
	switch (BioNet_edge.type) {
		case ACTIVATORY:
			out << " ACT ";
			break;
		case INHIBITORY:
			out << " REP ";
			break;
		case UNKNOWN:
			out << " UNKNOWN ";
			break;
	}
	return out;
}

Module::Module() {
	// TODO Auto-generated constructor stub
}

Module::~Module() {
	// TODO Auto-generated destructor stub
}

void Module::TestPrint(ostream* f) {
	cout << "{-- " << name << " --}" << endl;
	GeneCircuitGraph::TestPrint(f);
	*f << "Inputs: ";
	for (unsigned int i = 0; i < inputs.size(); i++)
		*f << inputs[i] << " ";
	*f << endl;
	*f << "Outputs: ";
	for (unsigned int i = 0; i < outputs.size(); i++)
		*f << outputs[i] << " ";
	*f << endl;
}

Module* GenerateGateNetwork(int number_of_gates, int number_of_inputs, bool is_complete) {
	ARGEdit* editor_tmp = new ARGEdit;
	int n1 = editor_tmp->InsertNode(new BioNetNode("Ligand",""));
	int n2 = editor_tmp->InsertNode(new BioNetNode("Protein","gfp"));
	editor_tmp->InsertEdge(n1,n2, new BioNetEdge("UNKNOWN"));
	vector<IdPair*> input_id_list; // each pair contains source and dest of the input edge
	input_id_list.push_back(new IdPair(n1,n2));

	for (int g = 0; g < number_of_gates - number_of_inputs; g++) {
		int id = rand() % input_id_list.size();
		int gate_type = rand() % 3;
		switch (gate_type) {
			case 0: {
				// NOT
				delete ((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[id]->first))->component;
				((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[id]->first))->component = CreateComponent("NOT_GATE");
				int new_id = editor_tmp->InsertNode(new BioNetNode("Ligand",""));
				editor_tmp->InsertEdge(new_id, input_id_list[id]->first, new BioNetEdge("UNKNOWN"));
				input_id_list[id]->second = input_id_list[id]->first;
				input_id_list[id]->first = new_id;
				break;
			}
			case 2: {
				// OR
				delete ((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[id]->first))->component;
				((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[id]->first))->component = CreateComponent("OR_GATE");
				int new_id_1 = editor_tmp->InsertNode(new BioNetNode("Ligand",""));
				int new_id_2 = editor_tmp->InsertNode(new BioNetNode("Ligand",""));
				editor_tmp->InsertEdge(new_id_1, input_id_list[id]->first, new BioNetEdge("UNKNOWN"));
				editor_tmp->InsertEdge(new_id_2, input_id_list[id]->first, new BioNetEdge("UNKNOWN"));
				input_id_list[id]->second = input_id_list[id]->first;
				input_id_list[id]->first = new_id_1;
				input_id_list.push_back(new IdPair(new_id_2, input_id_list[id]->second));
				break;
			}
			case 1: {
				// AND
				delete ((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[id]->first))->component;
				((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[id]->first))->component = CreateComponent("AND_GATE");
				int new_id_1 = editor_tmp->InsertNode(new BioNetNode("Ligand",""));
				int new_id_2 = editor_tmp->InsertNode(new BioNetNode("Ligand",""));
				editor_tmp->InsertEdge(new_id_1, input_id_list[id]->first, new BioNetEdge("UNKNOWN"));
				editor_tmp->InsertEdge(new_id_2, input_id_list[id]->first, new BioNetEdge("UNKNOWN"));
				input_id_list[id]->second = input_id_list[id]->first;
				input_id_list[id]->first = new_id_1;
				input_id_list.push_back(new IdPair(new_id_2, input_id_list[id]->second));
				break;
			}
		}
	}

	if (input_id_list.size() > (unsigned int)number_of_inputs) { // the number of inputs is larger than the expected one, so we need to merge them together
		vector<vector<int> > dest_id_matrix_tmp(number_of_inputs);	// store the dest node ids of each set of input nodes
		for (int id_tmp = 0; id_tmp < number_of_inputs; id_tmp++)
			dest_id_matrix_tmp[id_tmp].push_back(input_id_list[id_tmp]->second);
		vector<int> del_node_id_list;

		for (unsigned int id_gather = number_of_inputs; id_gather < input_id_list.size(); id_gather++) {
			int tmp;
			bool check;
			do {// check to ensure that there are no two input nodes that connect to the same dest node can be merged together
				tmp = rand() % number_of_inputs;
				check = true;
				for (unsigned int it = 0; it < dest_id_matrix_tmp[tmp].size(); it++)
					if (input_id_list[id_gather]->second == dest_id_matrix_tmp[tmp][it]) {
						check = false;
						break;
					}
			} while (!check);
			del_node_id_list.push_back(input_id_list[id_gather]->first);
			editor_tmp->InsertEdge(input_id_list[tmp]->first, input_id_list[id_gather]->second, new BioNetEdge("UNKNOWN"));
			dest_id_matrix_tmp[tmp].push_back(input_id_list[id_gather]->second);
			editor_tmp->DeleteEdge(input_id_list[id_gather]->first, input_id_list[id_gather]->second);
		}
	}

	for (unsigned int i = 0; i < input_id_list.size(); i++) { // Add an YES_GATE for each input node
		int id_tmp = editor_tmp->InsertNode(new BioNetNode("Ligand"));
		delete ((BioNetNode*)editor_tmp->GetNodeAttr(id_tmp))->component;
		((BioNetNode*)editor_tmp->GetNodeAttr(id_tmp))->component = ((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[i]->first))->component;
		((BioNetNode*)editor_tmp->GetNodeAttr(input_id_list[i]->first))->component = CreateComponent("YES_GATE");
		editor_tmp->InsertEdge(id_tmp,input_id_list[i]->first,new BioNetEdge("UNKNOWN"));
	}
	bool del_check;
	// remove abundant nodes
	do {
		del_check = true;
		for (int i = 0; i < editor_tmp->NodeCount(); i++)
			if (editor_tmp->OutEdgeCount(i) == 0 && ((BioNetNode*)editor_tmp->GetNodeAttr(i))->component->getType() != PROTEIN) {
				editor_tmp->DeleteNode(i);
				del_check = false;
				break;
		}
	}
	while (!del_check);
	return new Module(editor_tmp);
}
