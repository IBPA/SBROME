/*
 * DesignProblem.h
 *
 *  Created on: Apr 20, 2012
 *      Author: linh, UC Davis
 */

#ifndef DESIGNPROBLEM_H_
#define DESIGNPROBLEM_H_

#include <iostream>
#include <string>
#include "GeneCircuitGraph.h"


class SimulationData {
public:
	IdList inputs, outputs;						// List of ids of input nodes & output nodes
	ValueMatrix input_values, output_values;	// number of rows = number of data points (for steady state cases) or number of time points (temporal cases)
	void TestPrint(ostream* f);
};

class DesignProblem {
public:
	DesignProblem();
	DesignProblem(string filename);
	DesignProblem(istream* f);
	DesignProblem(const DesignProblem& another_design_problem);
	DesignProblem& operator = (const DesignProblem& another_design_problem);

	void TestPrint(ostream* f);
	GeneCircuitGraph gcg;
	SimulationData desired_behavior;

private:
	void Init(istream* f);
	//void ClearAll();
};

#endif /* DESIGNPROBLEM_H_ */
