#include <map>
//#include <iostream>
#include <fstream>
#include "json/json.h"


using namespace std;

void readConstraints(std::map<char, int> *constraintMap)
{	
	const char* filename = "constraints.json";

	//constraintMap->insert(std::pair<char,int>('a',100));
	//constraintMap->insert(std::pair<char,int>('b',140));

	Json::Value root;

	std::ifstream config_doc(filename, std::ifstream::binary);

	config_doc >> root;

	const Json::Value user_constraints = root["constraints"];

	for (int i=0; i<user_constraints.size(); i++){
		
		//printf("%s \n",user_constraints[i].asString().c_str());

		char c = user_constraints[i].asString().c_str()[0];

		constraintMap->insert(std::pair<char,int>(c,1));

	}
	

	return;


}