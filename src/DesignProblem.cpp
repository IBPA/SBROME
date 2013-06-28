/*
 * DesignProblem.cpp
 *
 *  Created on: Apr 20, 2012
 *      Author: linh, UC Davis
 */

#include <fstream>
#include "DesignProblem.h"

void SimulationData::TestPrint(ostream* f) {
	for (int row = 0; row < input_values.size(); row++) {
		for (int col = 0; col < input_values[row].size(); col++)
			*f << input_values[row][col] << "\t";
		for (int col = 0; col < output_values[row].size(); col++)
			*f << output_values[row][col] << "\t";
		*f << endl;
	}
}

DesignProblem::DesignProblem() {

}

DesignProblem::DesignProblem(string filename) {
	std::ifstream f;
	f.open(filename.c_str());
	if (!f.good()) {
		cout << "Can not open file: " << filename << endl;
		return;
	}
	Init(&f);
	f.close();
}

DesignProblem::DesignProblem(const DesignProblem& another_design_problem) {
	gcg = another_design_problem.gcg;
	desired_behavior = another_design_problem.desired_behavior;
}

DesignProblem::DesignProblem(istream* f) {
	Init(f);
}

DesignProblem& DesignProblem::operator = (const DesignProblem& another_design_problem) {
	gcg = another_design_problem.gcg;
	desired_behavior = another_design_problem.desired_behavior;
	return *this;
}

void DesignProblem::Init(istream *f) {
	/*
	gcg.Init(f);
	string comments;
	*f >> comments;
	int number_of_conditions;
	*f >> number_of_conditions;
	inputs.resize(number_of_conditions);
	outputs.resize(number_of_conditions);
	for (int cond = 0; cond < number_of_conditions; cond++) {
		inputs[cond].resize(gcg.getNumberOfInputs());
		outputs[cond].resize(gcg.getNumberOfOutputs());
		for (int i = 0; i < gcg.getNumberOfInputs(); i++)
			*f >> inputs[cond][i];
		for (int i = 0; i < gcg.getNumberOfOutputs(); i++)
			*f >> outputs[cond][i];
	}
	*/
}

void DesignProblem::TestPrint(ostream* f){
	gcg.TestPrint(f);
	*f << "INPUTS: " ;
	for (unsigned int i = 0; i < desired_behavior.inputs.size(); i++)
		*f << desired_behavior.inputs[i] << "\t";
	*f << endl << "OUTPUTS: ";
	for (unsigned int i = 0; i < desired_behavior.outputs.size(); i++)
		*f << desired_behavior.outputs[i] << "\t";
	*f << endl;
	for (unsigned int i = 0; i < desired_behavior.input_values.size(); i++) {
		for (int j = 0; j < desired_behavior.input_values[i].size(); j++)
			*f << desired_behavior.input_values[i][j] << "\t";
		for (int j = 0; j < desired_behavior.output_values[i].size(); j++)
			*f << desired_behavior.output_values[i][j] << "\t";
		*f << endl;
	}
}
