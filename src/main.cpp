// Simple main.cpp
// for perfect matching of bonds from XYZ files

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <stdio.h>

#include "tokenize.h"
#include "molecule.h"
#include "atom.h"
#include "bond.h"

bool readXYZ(Molecule &mol, const char* filename);
void readConstraints(std::map<char, int> *constraintMap);

using namespace std;

int main (int argc, char *argv[])
{
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " {filename}" << endl;
    return 1;
  }

  Molecule mol;
  // loop through the filenames
  for (unsigned int a = 1; a < argc; ++a) {
    if (!readXYZ(mol, argv[a]))
      cout << "Cannot read the XYZ file" << endl;

    cout << " Molecule has " << mol.numberOfAtoms() << " atoms." << endl;
  }

  // This is the map for the constraints.
  std::map<char, int> constraint_map;
  readConstraints(&constraint_map);

  //Printing out all the constraints that were read from the file.
  for(auto it= constraint_map.begin(); it != constraint_map.end(); ++it){
    printf(" %c : %d \n", it->first, it->second);

  }
  

  


  return 0;
}
