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
 double prime = 512009;
 double err = 0.000001;
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
 void inverse(std::vector<std::vector<double>> &C, int N, std::map<int, int> excl);
 double determinant(std::vector<std::vector<double>> A, int N);
 void printMatrix(std::vector<std::vector<double>> A);
 void printAugMatrix();
 void printDeterminant();
 std::vector<std::vector<double>> copy_graph();
 std::vector<std::vector<double>> copy_augC();
 std::vector<std::tuple<int, int>> matching();


 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;
  unsigned int numBonds;
  std::vector<std::vector<double>> graph;
  std::vector<std::vector<double>> augC;
};

#endif
