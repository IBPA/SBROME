/*
 * GRN.h
 *
 *  Created on: Apr 8, 2012
 *      Author: linh, UC Davis
 */

#ifndef GRN_H_
#define GRN_H_

#include <vector>
#include <sstream>
#include <string>
#include <map>

#ifndef NULL
	#define NULL	0
#endif

#define DEBUG

using namespace std;

typedef std::vector<int> IdList;
typedef std::vector<IdList> IdMatrix;
typedef std::pair<int,int> IdPair;
typedef std::vector<IdPair> IdPairList;
typedef std::vector<double> ValueList;
typedef std::vector<ValueList> ValueMatrix;
typedef std::vector<std::string> StringList;
typedef std::vector<StringList> StringMatrix;
typedef std::map<std::string, std::pair<StringList,int> > StrTree;
typedef std::map<std::string,int> CheckList;

enum LEVEL {
	HIGH,
	LOW
};

enum FinalStatus {
	SUCCESS,
	FAILURE
};

#endif /* GRN_H_ */
