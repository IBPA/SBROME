/*
 * GeneCircuitGraph.h
 *
 *  Created on: May 2, 2012
 *      Author: linh, UC Davis
 */

#ifndef GENECIRCUITGRAPH_H_
#define GENECIRCUITGRAPH_H_

#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "GRN.h"
#include "VFLIB/argraph.h"
#include "VFLIB/argloader.h"
#include "VFLIB/allocpool.h"

#define free_var(x) if (x != NULL) {delete x; x = NULL;}
#define free_array(x) if (x != NULL) {delete [] x; x = NULL;}

enum ReturnValue {
	STR_EQ	= 0
};

enum Category {
	CIRCUIT_SIGNAL,
	LOGIC_GATE,
	MOLECULE,
};

enum IOSignal{
	INTERMEDIATE = 1000,
	INPUT = 1001,
	OUTPUT = 1002,
	LINKER = 1003
};

enum GateType {
	UNKNOWN_GATE 		= -1,
	DUMMY_GATE			= 0,
	YES_GATE			= 100,
	NOT_GATE 			= 101,
	AND_GATE 			= 102,
	NAND_GATE 			= 103,
	OR_GATE 			= 104,
	NOR_GATE			= 105,
	OSCILLATOR			= 106

};

enum MoleculeType {
	UNKNOWN_COMPLEX_MOLECULE 	= -2,
	UNKNOWN_MOLECULE 			= -1,
	LIGAND 						= 0,
	m_RNA 						= 1,
	t_RNA						= 2,
	PROTEIN 					= 3,
	LIGAND_COMPLEX 				= 4,
	RNA_COMPLEX 				= 5,
	PROTEIN_COMPLEX 			= 6,
	LIGAND_PROTEIN_COMPLEX 		= 7,
	POOL 						= 8
};

enum mRNASubType {
	REGULAR	= 1,
	TANDEM 	= 2,
	HYBRID	= 3
};

enum ProteinSubType {
	REPORTER	= 1
};

enum InteractionType {
	UNKNOWN 		= -1,
	ACTIVATORY		= 0,
	INHIBITORY		= 1,
	ACTIVATORY_OR	= 2,
	INHIBITORY_OR	= 3,
	NONE			= 4
};

string Id2Str(int type_id);
int Str2Id(string type_str);
typedef std::pair<int,int> IdPair;

struct NodePair {
	NodePair(node_id n1, node_id n2) {
		small_graph_node = n1;
		large_graph_node = n2;
	}
	node_id small_graph_node;
	node_id large_graph_node;
};

struct Link {
	Link (int src_, int dest_) {
		src = src_;
		dest = dest_;
	}
	int src;
	int dest;
};

typedef vector<NodePair> GraphIsomorphism;
typedef vector<GraphIsomorphism> GraphIsomorphismCollection;
typedef vector<Link> ListOfLinks;

class GeneCircuitComponent {
public:
	GeneCircuitComponent() {
		name = "";
		variant_id = UNKNOWN;
	}
	virtual ~GeneCircuitComponent(){
		// Do nothing
	}
	virtual bool CompareTo(GeneCircuitComponent* another_component) = 0; // circuit.CompareTo(module)
	virtual int getCategory() = 0;
	virtual int getType()= 0;
	virtual void operator >> (istream& in) = 0;
	virtual void operator << (ostream& out) = 0;
	virtual GeneCircuitComponent* clone() = 0;

	string getName() {
		return name;
	}

	string name;
	int variant_id;
};

class Molecule: public GeneCircuitComponent {// Each node represents a molecular species or a compounds of some molecular species in a biological network
public:
	int getCategory() {
		return MOLECULE;
	}
	int getType() {
		return type;
	}

protected:
	int type;	// Ligand, protein, ligand-protein, etc
};

class SingleMolecule: public Molecule {
public:
	SingleMolecule() {
		type = UNKNOWN_MOLECULE;
		sub_type = UNKNOWN;
		name = "";
	}
	SingleMolecule(int type_, string name_, int sub_type_ = UNKNOWN) {
		type = type_;
		name = name_;
		sub_type = sub_type_;
	}
	virtual ~SingleMolecule() {
		// Do nothing
	}

	bool CompareTo(GeneCircuitComponent* another_component);

	void operator >> (istream& in) {
		// Type is already determined before that
		// Name
		in >> name;
		if (name.compare("UNKNOWN") == STR_EQ)
			name = "";
	}

	void operator << (ostream& out) {
		out << Id2Str(type) << "\t" << name;
	}

	GeneCircuitComponent* clone() {
		return new SingleMolecule(type, name);
	}
public:
	int sub_type;	// mRNA: hybrid, tandem; Protein: reporter
};

class ComplexMolecule: public Molecule {	// complex of two single molecules
public:
	ComplexMolecule() {
		first_mol = second_mol = NULL;
		type = UNKNOWN_COMPLEX_MOLECULE;
		name = UNKNOWN;
	}
	ComplexMolecule(const ComplexMolecule &another_molecule) {
		type = another_molecule.type;
		name = another_molecule.name;
		first_mol = (SingleMolecule*)another_molecule.first_mol->clone();
		second_mol = (SingleMolecule*)another_molecule.second_mol->clone();
	}
	ComplexMolecule(int type_, string name_) {
		type = type_;
		name = name_;
		switch (type) {
			case LIGAND_PROTEIN_COMPLEX:
				first_mol = new SingleMolecule(LIGAND, "");
				second_mol = new SingleMolecule(PROTEIN, "");
				break;
			case RNA_COMPLEX:
				first_mol = new SingleMolecule(m_RNA, "");
				second_mol = new SingleMolecule(m_RNA, "");
				break;
			case PROTEIN_COMPLEX:
				first_mol = new SingleMolecule(PROTEIN, "");
				second_mol = new SingleMolecule(PROTEIN, "");
				break;
		}
	}
	ComplexMolecule(SingleMolecule* m1, SingleMolecule* m2) {
		// First
		if (m1 != NULL)
			first_mol = (SingleMolecule*)m1->clone();
		else
			first_mol = NULL;
		// Second
		if (m2 != NULL)
			second_mol = (SingleMolecule*)m2->clone();
		else
			second_mol = NULL;
		// Type and Name
		if (m1 != NULL && m2 != NULL) {
			if (m1->getType() == LIGAND)
				type = LIGAND_PROTEIN_COMPLEX;
			else if (m1->getType() == m_RNA)
				type = RNA_COMPLEX;
			else if (m1->getType() == PROTEIN)
				type = PROTEIN_COMPLEX;
			name = m1->getName() + "-" + m2->getName();
		}
		else {
			type = UNKNOWN_COMPLEX_MOLECULE;
			name = "";
		}
	}
	virtual ~ComplexMolecule() {
		if (first_mol != NULL)
			delete first_mol;
		if (second_mol != NULL)
			delete second_mol;
	}
	bool CompareTo(GeneCircuitComponent* another_component) { // circuit.CompareTo(module)
		if (another_component->getType() != this->getType())
			return false;
		else
			return true;
	}
	void operator >> (istream& in) {
		// Type
		string str_type;
		in >> str_type;
		if (str_type.compare("LIGAND_PROTEIN_COMPLEX") == STR_EQ) {
			type = LIGAND;
		}
		else if (str_type.compare("RNA_COMPLEX") == STR_EQ) {
			type = m_RNA;
		}
		else if (str_type.compare("PROTEIN_COMPLEX") == STR_EQ) {
			type = PROTEIN;
		}
		first_mol = new SingleMolecule();
		second_mol = new SingleMolecule();
		first_mol->operator >>(in);
		second_mol->operator >>(in);
		name = first_mol->getName().append(second_mol->getName());
	}
	void operator << (ostream& out) {
		if (name.compare("UNKNOWN") != STR_EQ)
			out << Id2Str(type) << "\t" << name;
		else
			out << Id2Str(type) << "\t" << first_mol->getName() << "-" << second_mol->getName();
	}

	GeneCircuitComponent* clone() {
		return new ComplexMolecule(*this);
	}

private:
	SingleMolecule *first_mol, *second_mol;
};

class Device: public GeneCircuitComponent {
	string motif_name;				// use for organizing and visualizing
	vector<string> node_name_list;	// this list is sorted by the node indexes of the motif above
};

class OneOutputLogicGate: public Device {
public:
	int getCategory() {
		return LOGIC_GATE;
	}
};

class DummyGate: public OneOutputLogicGate {	// Only for the scalability check
public:
	DummyGate() {
		gate_type = UNKNOWN;
		gate_topology = UNKNOWN;
	}
	DummyGate(int gate_type_) {
		gate_type = gate_type_;
		gate_topology = UNKNOWN;
	}
	DummyGate(int gate_type_, int gate_topology_, string name_) {
		gate_type = gate_type_;
		gate_topology = gate_topology_;
		name = name_;
	}
	DummyGate(const DummyGate &another_gate) {
		gate_type = another_gate.gate_type;
		gate_topology = another_gate.gate_topology;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << Id2Str(gate_type) << "\t" << name;
	}
	int getType() {
		return DUMMY_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) {
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == DUMMY_GATE && (gate_type == ((DummyGate*)another_component)->gate_type)
			&& (gate_topology == UNKNOWN || gate_topology == ((DummyGate*)another_component)->gate_topology)
			&& (name.compare("") == STR_EQ || name.compare(another_component->name) == STR_EQ))
				return true;
		else
			return false;
	}
	GeneCircuitComponent* clone() {
		return new DummyGate(*this);
	}

	// name is inherited from the GeneCircuitComponent class
	int gate_type;	// YES_GATE, NOT_GATE, OR_GATE or AND_GATE
	int gate_topology;
};

class OneInputGate: public OneOutputLogicGate {
protected:

};

class YESGate: public OneInputGate {
public:
	YESGate() {
		// Do nothing
	}
	YESGate(const YESGate &another_gate) {
		name = another_gate.name;
	}
	virtual ~YESGate() {
		// Do nothing
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "YESGate \t" << name;
	}
	int getType() {
		return YES_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) { // circuit.CompareTo(module)
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == YES_GATE
			&& (name.compare("") == STR_EQ || name.compare(another_component->name) == STR_EQ))
				return true;
		else
			return false;
	}
	GeneCircuitComponent* clone() {
		return new YESGate(*this);
	}
};

class NOTGate: public OneInputGate {
public:
	NOTGate() {
		// Do nothing
	}
	NOTGate(const NOTGate &another_gate) {
		name = another_gate.name;
	}
	virtual ~NOTGate() {
		// Do nothing
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "NOTGate \t" << name;
	}
	int getType() {
		return NOT_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) { // circuit.CompareTo(module)
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == NOT_GATE
		&& (name.compare("") == STR_EQ || name.compare(another_component->name) == STR_EQ))
			return true;
		else
			return false;
	}
	GeneCircuitComponent* clone() {
		return new NOTGate(*this);
	}
};

class TwoInputGate: public OneOutputLogicGate {
protected:
};

class ANDGate: public TwoInputGate {
public:
	ANDGate() {
		// Do nothing
	}
	ANDGate(const ANDGate &another_gate) {
		name = another_gate.name;
	}
	virtual ~ANDGate() {
		// Do nothing
	}
	int getType() {
		return AND_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) { // circuit.CompareTo(module)
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == AND_GATE
		&& (name.compare("") == STR_EQ || name.compare(another_component->name)))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "ANDGate \t" << name;
	}
	GeneCircuitComponent* clone() {
		return new ANDGate(*this);
	}
};

class NANDGate: public TwoInputGate {
public:
	NANDGate() {
		// Do nothing
	}
	NANDGate(const NANDGate &another_gate) {
		name = another_gate.name;
	}
	virtual ~NANDGate() {
		// Do nothing
	}
	int getType() {
		return NAND_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) {
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == NAND_GATE
		&& (name.compare("") == STR_EQ || name.compare(another_component->name)))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "NANDGate \t" << name;
	}
	GeneCircuitComponent* clone() {
		return new NANDGate();
	}
};

class ORGate: public TwoInputGate {
public:
	ORGate() {
		// Do nothing
	}
	ORGate(const ORGate &another_gate) {
		name = another_gate.name;
	}
	virtual ~ORGate() {
		// Do nothing
	}
	int getType() {
		return OR_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) {
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == OR_GATE
		&& (name.compare("") == STR_EQ || name.compare(another_component->name)))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "ORGate \t" << name;
	}
	GeneCircuitComponent* clone() {
		return new ORGate(*this);
	}
};

class NORGate: public TwoInputGate {
public:
	NORGate() {
		// Do nothing
	}
	NORGate(const NOTGate &another_gate) {
		name = another_gate.name;
	}
	virtual ~NORGate() {
		//ClearAll();
	}
	int getType() {
		return NOR_GATE;
	}
	bool CompareTo(GeneCircuitComponent* another_component) {
		if (another_component->getCategory() == LOGIC_GATE && another_component->getType() == NOR_GATE
		&& (name.compare("") == STR_EQ || name.compare(another_component->name)))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "NORGate \t" << name;
	}
	GeneCircuitComponent* clone() {
		return new NORGate(*this);
	}
};

class CircuitSignal: public GeneCircuitComponent {
public:
	int getCategory() {
		return CIRCUIT_SIGNAL;
	}
	GeneCircuitComponent* physical_instance;
};

class InputSignal: public CircuitSignal {
public:
	InputSignal() {
		physical_instance = NULL;
	}
	InputSignal(const InputSignal& another_input) {
		name= another_input.name;
		if (another_input.physical_instance != NULL)
			physical_instance = another_input.physical_instance->clone();
		else
			physical_instance = NULL;
	}
	~InputSignal() {
		free_var(physical_instance);
	}
	int getType() {
		return INPUT;
	}
	bool CompareTo(GeneCircuitComponent* another_component) {
		if (another_component->getCategory() == CIRCUIT_SIGNAL && another_component->getType() == INPUT && (physical_instance == NULL || physical_instance->CompareTo(((InputSignal*) another_component)->physical_instance)))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "INPUT \t" << name << "\t";
		//if (physical_instance != NULL)
		//	physical_instance->operator <<(out);
	}

	GeneCircuitComponent* clone() {
		return new InputSignal(*this);
	}
};

class OutputSignal: public CircuitSignal {
public:
	OutputSignal() {
		physical_instance = NULL;
	}
	OutputSignal(const OutputSignal& another_output) {
		name= another_output.name;
		if (another_output.physical_instance != NULL)
			physical_instance = another_output.physical_instance->clone();
		else
			physical_instance = NULL;
	}
	~OutputSignal() {
		free_var(physical_instance);
	}
	int getType() {
		return OUTPUT;
	}
	bool CompareTo(GeneCircuitComponent* another_component) {
		if (another_component->getCategory() == CIRCUIT_SIGNAL && another_component->getType() == OUTPUT
		    && (physical_instance == NULL || physical_instance->CompareTo(((InputSignal*) another_component)->physical_instance)))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "OUTPUT \t" << name << "\t";
		//if (physical_instance != NULL)
		//	physical_instance->operator <<(out);
	}
	GeneCircuitComponent* clone() {
		return new OutputSignal(*this);
	}
};

class LinkSignal: public CircuitSignal {
public:
	LinkSignal() {
		physical_instance = NULL;
	}
	LinkSignal(const LinkSignal& another_linker) {
		name = another_linker.name;
		if (another_linker.physical_instance != NULL)
			physical_instance = another_linker.physical_instance->clone();
		else
			physical_instance = NULL;
	}
	~LinkSignal() {
		free_var(physical_instance);
	}
	int getType() {
		return LINKER;
	}
	bool CompareTo(GeneCircuitComponent* another_component) { // circuit.CompareTo(module)
		if ((another_component->getCategory() == CIRCUIT_SIGNAL && (physical_instance == NULL || physical_instance->CompareTo(((InputSignal*) another_component)->physical_instance)))
			|| (another_component->getCategory() == MOLECULE && (physical_instance == NULL || physical_instance->CompareTo(another_component))))
			return true;
		else
			return false;
	}
	void operator >> (istream& in) {

	}
	void operator << (ostream& out) {
		out << "LINKER \t" << name << "\t";
		//if (physical_instance != NULL)
		//	physical_instance->operator <<(out);
	}
	GeneCircuitComponent* clone() {
		return new LinkSignal(*this);
	}
};

GeneCircuitComponent* CreateComponent (string type);
GeneCircuitComponent* CreateComponent (string type, string name);

// BioNetNode: purely abstract class
class BioNetNode {
public:
	GeneCircuitComponent* component;

	~BioNetNode() {
		free_var(component);
	}
	BioNetNode() {
		component = NULL;
		original_node_id = 0;
	}
	BioNetNode(string type) {
		component = CreateComponent(type);
		component->name = "";
		original_node_id = 0;
	}
	BioNetNode(string type, string name) {
		component = CreateComponent(type);
		component->name = name;
		original_node_id = 0;
	}
	BioNetNode(const BioNetNode& another_node) {
		component = another_node.component->clone();
		module_id_list = another_node.module_id_list;
		original_node_id = another_node.original_node_id;
	}
	BioNetNode* clone() {
		return new BioNetNode(*this);
	}

	vector<int> module_id_list;	// it is used for the module matching step
	int original_node_id;
};

BioNetNode* CreateBioNetNode (string type);
BioNetNode* CreateBioNetNode (string type, string name);

class BioNetNodeDestroyer: public AttrDestroyer {
public:
	virtual void destroy (void *p) {
		delete ((BioNetNode*)p);
	}
	virtual BioNetNodeDestroyer* clone() {
		return new BioNetNodeDestroyer();
	}
};

class BioNetNodeComparator: public AttrComparator {
	virtual bool compatible(void* pa, void* pb) {	// pa: module, pb: circuit
		BioNetNode* module = (BioNetNode*) pa;
		BioNetNode* circuit = (BioNetNode*) pb;
		if (!circuit->module_id_list.empty() && (module->component->getCategory() != CIRCUIT_SIGNAL || module->component->getType() != INPUT))
			return false;
		else
			return ((BioNetNode*)pb)->component->CompareTo(((BioNetNode*)pa)->component);
	}
	virtual BioNetNodeComparator* clone() {
		return new BioNetNodeComparator();
	}
};

class BioNetEdge {
public:
	BioNetEdge() {
		// Do nothing
		type = UNKNOWN;
	}
	BioNetEdge(int type_) {
		type = type_;
	}
	BioNetEdge(string relation_type) {
		if (relation_type.compare("ACT") == STR_EQ || relation_type.compare("ACT_OR") == STR_EQ)
			type = ACTIVATORY;
		else if (relation_type.compare("REP") == STR_EQ || relation_type.compare("REP_OR") == STR_EQ)
			type = INHIBITORY;
		else if (relation_type.compare("UNKNOWN") == STR_EQ)
			type = UNKNOWN;
		else
			cout << "Relationship: " << relation_type << " has not defined yet";
	}

	BioNetEdge(BioNetEdge* another_edge) {
		type = another_edge->type;
	}
	BioNetEdge* clone() {
		return new BioNetEdge(this);
	}
	int type;
	//int src_molecule_node_id;		// id of the molecule in the source GC module
	//int dest_molecule_node_id;		// id of the molecule in the destination GC module

};

class BioNetEdgeDestroyer: public AttrDestroyer {
public:
	virtual void destroy (void *p) {
		delete ((BioNetEdge*)p);
	}
	virtual BioNetEdgeDestroyer* clone() {
		return new BioNetEdgeDestroyer();
	}
};

class BioNetEdgeComparator: public AttrComparator {
	virtual bool compatible(void* pa, void* pb) {
		BioNetEdge* module_edge = (BioNetEdge*) pa;	// module
		BioNetEdge* circuit_edge = (BioNetEdge*) pb;	// circuit
		return ((circuit_edge->type == UNKNOWN) || (module_edge->type == circuit_edge->type));
	}
	virtual BioNetEdgeComparator* clone() {
		return new BioNetEdgeComparator();
	}
};

// Forward declaration
class PartDatabase;

class GeneCircuitGraph: public ARGraph<BioNetNode, BioNetEdge> {
public:
	GeneCircuitGraph();
	GeneCircuitGraph(ARGLoader* loader);
	GeneCircuitGraph(istream* f);
	GeneCircuitGraph(const GeneCircuitGraph& another_gene_circuit_graph);
	GeneCircuitGraph& operator = (const GeneCircuitGraph& another_gene_circuit_graph);
	virtual ~GeneCircuitGraph();

	ARGLoader* ConvertToLoader();
	void TestPrint(ostream* f);
	void PrintAll(PartDatabase* partDB, ostream* f);
	bool IsComplete();
	int getNumberOfInputs();
	int getNumberOfOutputs();
	int getNumberOfModuleNodes();
	int getNumberOfUnkownNodes();	// if all node is determined completely

	void Init(istream* f);
	GeneCircuitGraph* clone();
private:
	int number_of_inputs, number_of_outputs;
};

class Module: public GeneCircuitGraph {
public:
	Module();
	Module(ARGLoader* loader): GeneCircuitGraph(loader){
		for (int i = 0; i < n; i++) {
			if (this->InEdgeCount(i) == 0)
				inputs.push_back(i);
			if (this->OutEdgeCount(i) == 0)
				outputs.push_back(i);
		}
	}
	Module(ARGLoader* loader, IdList inputs_, IdList outputs_): GeneCircuitGraph(loader){
		inputs = inputs_;
		outputs = outputs_;
	}
	virtual ~Module();
	void TestPrint(ostream* f);
	vector<string> getOutputName() {
		vector<string> tmp;
		for (unsigned int i = 0; i < outputs.size();i++)
			tmp.push_back(this->GetNodeAttr(outputs[i])->component->name);
		return tmp;
	}
	//void MatchingAtOneNode(Module* input_module, int node_id, GraphIsomorphismCollection* matching_list);	// match the input_module so that it's output can match well with the node at node_id

	string name;
	IdList inputs, outputs;
};

struct EdgeInfo{
	int src,dest;
	BioNetEdge* edge;
};

class NodeVisualInfo {
public:
	NodeVisualInfo() {
		pos_x = 0;
		pos_y = 0;
		width = 0;
		height = 0;
		line_width = 1;
	}
	NodeVisualInfo(double x, double y, double w, double h, double l_w = 1) {
		pos_x = x;
		pos_y = y;
		width = w;
		height = h;
		line_width = l_w;
	}
	double pos_x, pos_y;
	double width, height;
	double line_width;
};

typedef map<int, NodeVisualInfo*> NodeVisualInfoMap;

Module* GenerateGateNetwork(int number_of_gates, int number_of_inputs, bool is_complete);	// is_complete: if there are ligands and YES_GATE or not

#endif /* GENECIRCUITGRAPH_H_ */
