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
#include <ctime>


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


    clock_t begin = clock();


    mol.initializeGraph();


    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    /////

    printf("\n\nTime taken by initialize graph = %f\n\n", elapsed_secs);


    begin = clock();


    mol.perceiveBonds();

    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    /////

    printf("\n\nTime taken by perceive bonds= %f\n\n", elapsed_secs);

    cout << "Molecule has " << mol.numberOfAtoms() << " atoms and " << mol.numberOfBonds() << " bonds." << endl;

    begin = clock();
    std::vector<std::tuple<int,int>> M = mol.CUDAMatching();
    mol.printMatching(M);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    /////

    printf("\n\nTime taken by matching function = %f\n\n", elapsed_secs);

    mol.addMatchedBonds(M);
    mol.checkMatching(M);

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
  /*for(auto it= constraint_map.begin(); it != constraint_map.end(); ++it){
    printf(" %c : %d \n", it->first, it->second);

  }*/

  return 0;
}
