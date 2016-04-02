// Simple main.cpp
// for perfect matching of bonds from XYZ files

#include <vector>
#include <string>

#include <iostream>

#include "tokenize.h"
#include "molecule.h"
#include "atom.h"
#include "bond.h"

bool readXYZ(Molecule &mol, const char* filename);

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

    mol.perceiveBonds();

    cout << " Molecule has " << mol.numberOfAtoms() << " atoms and " << mol.numberOfBonds() << " bonds." << endl;
  }

  return 0;
}
