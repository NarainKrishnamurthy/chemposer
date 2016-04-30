// Simple Molecule class

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <cmath> 

#include "atom.h"
#include "bond.h"
#include <map>

class Molecule {
 public:
 Molecule() {};
 ~Molecule() {}

 void addAtom(Atom *atom) { _atoms.push_back(atom); }
 void addAtom(unsigned int element, double x, double y, double z)
 { Atom *a = new Atom(element, x, y, z, _atoms.size()); _atoms.push_back(a); }
 void addBond(Bond *bond) { _bonds.push_back(bond); }
 void addBond(Atom *a, Atom *b, unsigned short type = 0)
 { Bond *bond = new Bond(a, b); bond->setType(type); _bonds.push_back(bond); }

 unsigned int numberOfAtoms() { return _atoms.size(); }
 unsigned int numberOfBonds() { return numBonds; }

 std::vector<Atom*> atoms() { return _atoms; }
 std::vector<Bond*> bonds() { return _bonds; }

 void clear() { _atoms.clear(); _bonds.clear(); }
 
 void setNumBonds(unsigned num){
 	numBonds = num;
 }

 void perceiveBonds();
 void doMatching();
 void printMolecule();
 void printGraph();
 void initializeGraph();
 void inverse(std::map<int, int> *a, std::map<int, int> *b, float *err);
 void determinant(float *err);
// void determinant(int N);
 void printAugMatrix();

 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;
  unsigned int numBonds;
  std::vector<std::vector<float>> graph;
  std::vector<std::vector<float>> augC;
  float det;

};

#endif
