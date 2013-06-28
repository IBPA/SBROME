/*
 * PartDatabase.h
 *
 *  Created on: Aug 8, 2012
 *      Author: linh
 */

#ifndef PARTDATABASE_H_
#define PARTDATABASE_H_

#include <string>
#include <sqlite3.h>
#include <sstream>
#include "VFLIB/argedit.h"
#include "GRN.h"
#include "GeneCircuitGraph.h"

struct QueryResultInfo {
	// For access DB
	int rc;
	char *zErrMsg;
	std::string query;
	StringMatrix string_matrix;
};

struct TranscriptionRegulation {
	int regulator_type;
	double binding_affinity;
	double cooperativity;
};

struct LigandRegulation {
	double dissociation;
	double cooperativity;
};

/*
class MoleculeInteraction {
public:
	MoleculeInteraction() {
		first_mol = NULL;
		second_mol = NULL;
	}
	~MoleculeInteraction() {
		delete first_mol;
		delete second_mol;
	}
	Molecule* first_mol;
	Molecule* second_mol;
	int interaction_type;
};
*/

class PartDatabase {
public:
	PartDatabase();
	PartDatabase(std::string db_filename);
	virtual ~PartDatabase();

	// Interface
	int getNumberOfProteins() const;
	BioNetNode* getBioNetNode(int module_id, int node_id) const;

	//double getDegradationRate(int protein_id) const;
	double getDegradationRate(string protein_name) const;
	//int getNumberOfMutants(int promoter_id) const;
	int getNumberOfMutants(string promoter_name) const;
	//std::string getProteinName(int protein_id) const;
	//std::string getPromoterName(int promoter_id) const;
	//std::string getInducerName(int inducer_id) const;
	//std::string getPromoterMutantName(int promoter_index, int mutant_index) const;

	//double getPromoterBasal(int promoter_index, int mutant_index) const;
	double getPromoterBasal(string promoter_name, int mutant_index) const;
	//double getPromoterStrength(int promoter_index, int mutant_index) const;
	double getPromoterStrength(string promoter_name, int mutant_index) const;
	//double getPromoterAffinity(int protein_index, int promoter_index) const;
	TranscriptionRegulation getPromoterRegulation(string regulator_name, string promoter_name) const;
	//double getPromoterCooperativity(int protein_index, int promoter_index) const;
	//int getTranscriptionFactorType(int protein_index, int promoter_index) const;

	//double getInducerDissociation(int protein_index, int inducer_index) const;
	//double getInducerCooperativity(int protein_index, int inducer_index) const;
	LigandRegulation getLigandRegulation(string inducer_name, string protein_name) const;

	vector<Module*>* RandomizeModuleLibrary(int library_size, int max_module_size);

	// For module library, load each module and integrate the inputs
	vector<Module*>* loadModuleLibrary();
	vector<Module*>* GenerateModuleLibrary(string filename);
	string getSpeciesType(string species_name);
	string getSpeciesInteractionType(string src_species_name, string dest_species_name);
	StringMatrix* getAllInteractionSpecies(string species_type);
	// For regulation
	int getRegulationType(Molecule* m1, Molecule* m2);
	// For functional module substitutions
	int getNumberOfSubsitution(string functional_module_type) const;
	Module* loadSubstitution(string functional_module_type, int derive_id) const;
	// For macro module substitutions
	Module* loadContent(string macro_module_type, string macro_module_name, vector<NodeVisualInfo> *module_visual_info_list = NULL) const;

private:
	sqlite3 *db;
	int number_of_proteins;
};

void Reset(QueryResultInfo &qri);
vector<string>* split_string(string str, char c);
void AppendStringMatrix(StringMatrix* sm1, StringMatrix* sm2);
string ReadSentence(ifstream* f);
QueryResultInfo Init();
string Num2StrPlus1(int num);
string Num2Str(int num);
string Num2Str(double num);

#endif /* PARTDATABASE_H_ */
