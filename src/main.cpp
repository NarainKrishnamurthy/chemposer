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
bool writeSDF(Molecule &mol, const char* filename);
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

    mol.initializeGraph();
  
    mol.perceiveBonds();

    printf(" \n\n Augmented Matrix before taking the inverse \n\n");

    mol.printAugMatrix();

    cout << "Molecule has " << mol.numberOfAtoms() << " atoms and " << mol.numberOfBonds() << " bonds." << endl;

    std::map<int, int> excl_cols;

    std::map<int, int> excl_rows;
/*
    for(int i=0; i<mol.numberOfAtoms(); i++){
      excl_rows.insert(std::pair<int,int>(i,1));
      excl_cols.insert(std::pair<int,int>(i,1));
  
    }
*/
    mol.inverse(&excl_cols, &excl_rows);

    printf(" \n\n Augmented Matrix after taking inverse \n\n");

    mol.printAugMatrix();

    //mol.printGraph();
    mol.doMatching();

    std::vector<Atom*> atoms = mol.atoms();
    unsigned int j = 0;
    for (std::vector<Atom*>::iterator i = atoms.begin(); i < atoms.end(); ++i, ++j) {
      if ((*i)->atomicNum() == 6 && (*i)->numberOfDoubleBonds() != 1) {
        cout << " failed matching on atom " << j << endl;
        break;
      }
    } // end check loop

    // write an SD output
    char *filename = argv[a];
    // change extension
    char *pExt = strrchr(filename, '.');
    if (pExt != NULL)
      strcpy(pExt, ".sdf");
    else
      strcat(filename, ".sdf");

    if (!writeSDF(mol, filename))
      cout << "Cannot write the SDF file" << endl;

  } // end (loop through command-line args)

  // This is the map for the constraints.
  std::map<char, int> constraint_map;
  readConstraints(&constraint_map);

  //Printing out all the constraints that were read from the file.
  for(auto it= constraint_map.begin(); it != constraint_map.end(); ++it){
    printf(" %c : %d \n", it->first, it->second);

  }

  return 0;
}
