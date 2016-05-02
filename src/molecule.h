// Simple Molecule class

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <cmath> 

#include "atom.h"
#include "bond.h"
#include <map>
#include "Eigen/Dense"

class Molecule {
 public:
 double prime = 32768011;
 double err = 0.000000001;
 //double err = 0.0000001;
 Molecule() {};
 ~Molecule() {}

 void printMatching(std::vector<std::tuple<int, int>> M);
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
 void printMatrix(std::vector<std::vector<double>> A);
 void printAugMatrix();
 void printDeterminant();
 std::vector<std::tuple<int, int>> matching();


 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;
  unsigned int numBonds;
  std::vector<std::vector<double>> graph;
  std::vector<std::vector<double>> augC;
};

#endif
