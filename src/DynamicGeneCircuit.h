/*
 * DynamicGeneCircuit.h
 *
 *  Created on: May 2, 2012
 *      Author: linh
 */

#ifndef DYNAMICGENECIRCUIT_H_
#define DYNAMICGENECIRCUIT_H_

#include <vector>
#include "PartDatabase.h"
#include "GeneCircuitGraph.h"
#include "GA/GASimpleGA.h"			// we're going to use the simple GA
#include "GA/GA1DArrayGenome.h" 	// and the 2D binary string genome
#include "GA/GABin2DecGenome.h"

#define RESOLUTION 8

enum SimulationStatus {
	CONVERGENCE,
	OVER_TIME
};
/*
enum SignalType {
	INPUT,
	OUTPUT,
	LINKAGE,
	INSIDE
};
*/
struct SimulationResult {
	SimulationStatus status;
	double error;
};

struct Bound {
	Bound (double l, double u);
	double lower;
	double upper;
};

// Forward declaration
class MoleculePool;

class Component {
public:
	Component();
	virtual ~Component();

	vector<Component*> sub_component_list;
	vector<MoleculePool*> linkages;						// the source node of all linkages that connect sub-component together
														// EMPTY iff it is not a complex component
	vector<MoleculePool*> inputs;						// Outside the module, except the global inputs
	vector<MoleculePool*> outputs;						// Inside the module

	virtual void operator << (ostream& out);
	void Linearize(vector<MoleculePool*>& pool_list);
	int variant_id;
	int sub_circuit_id;
private:
	virtual void ClearAll();
};

class MoleculePool: public Component {
public:
	MoleculePool();
	virtual ~MoleculePool();
	double SteadyStateChange(const PartDatabase* pdb);
	double EstimateSteadyStateBound(const PartDatabase* pdb);

	Molecule* molecule;
	vector<MoleculePool*> parents;
	vector<MoleculePool*> children;
	ValueList min_concentration_list, max_concentration_list, best_concentration_list;
	double concentration;
	double min_concentration, max_concentration;				// Only the temporary values
	ValueList current_interval_index;							// Should be float instead of integer
	void operator << (ostream& out);
};

class DynamicGeneCircuit {
public:
	DynamicGeneCircuit();
	DynamicGeneCircuit(Module* gene_circuit_graph);
	DynamicGeneCircuit(Module* gene_circuit_graph, int number_of_sub_circuits);
	virtual ~DynamicGeneCircuit();

	void PrintOut(ostream* f);
	SimulationStatus SimulateSteadyState(const PartDatabase* pdb, const ValueMatrix& input_values, ValueMatrix& output_values, double threshold_value = 1e-7, int max_iteration = 10000);
	SimulationStatus PartialSimulateSteadyState(const PartDatabase* pdb, int sub_circuit_id, double threshold_value = 1e-7, int max_iteration = 10000);
	SimulationStatus EstimateSteadyStateBound(const PartDatabase* pdb, const ValueMatrix& input_values, ValueMatrix& linkage_lb, ValueMatrix& linkage_ub, double threshold_value = 1e-7, int max_iteration = 10000);
	void GASearch(const PartDatabase* pdb, const ValueMatrix& input_values, const ValueMatrix& output_values);
	void ExhaustiveSearch(const PartDatabase* pdb, const ValueMatrix& input_values, const ValueMatrix& output_values);
	void ModularSearch(const PartDatabase* pdb, const ValueMatrix& input_values, const ValueMatrix& output_values);

	Component* root_component;
private:
	vector<MoleculePool*> molecule_pool_list;
	void RecursiveClear(Component* component_node);
};

class DummyForGA{										// Dummy class for GA search
public:
	static float GAFitness(GAGenome &g);				// Estimate the fitness of each promoter mutant combination for the GA search procedure
	static DynamicGeneCircuit* dynamic_gene_circuit;
	static const PartDatabase* part_database;
	//static const DesignProblem* design_problem;
	static const ValueMatrix* inputs;
	static const ValueMatrix* outputs;
	static vector<MoleculePool*> promoter_list;
	static IdList number_of_mutants_list;
};

MoleculePool* BioNetNode2MoleculePool(BioNetNode* node_tmp);
#endif /* DYNAMICGENECIRCUIT_H_ */
