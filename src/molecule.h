// Simple Molecule class

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <cmath>

#include "atom.h"
#include "bond.h"
//#include "match.h"
#include <map>
#include "Eigen/Dense"

using namespace Eigen;

class Molecule {
 public:
 double prime = 32768011;
 double err = 0.0000000001;
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
 unsigned int numberOfBonds() { return _bonds.size(); }

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
 std::vector<std::tuple<int, int>> CUDAMatching();
 void printMatrix(std::vector<std::vector<double>> A);
 void printAugMatrix();
 void printDeterminant();
 std::vector<std::tuple<int, int>> matching();
 VectorXd solve(MatrixXd A,  VectorXd b, int m);
 void addMatchedBonds(std::vector<std::tuple<int, int>> M);

 int getfirst(int i){
     Bond *b = _bonds[i];
     int firstID = b->firstID();
     return firstID;
 }


 int getsecond(int i){
     Bond *b = _bonds[i];
     int secondID = b->secondID();
     return secondID;
 }

 void update(int i){
    _bonds[i]->setType(2);
 }

 void printBonds(){
     for(auto it= _bonds.begin(); it != _bonds.end(); ++it){
         printf("%d \n", (*it)->getType());
     }
 }


 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;
  double *cudaGraph;
  unsigned int numBonds;
  std::vector<std::vector<double>> graph;
  std::vector<std::vector<double>> augC;
};

#endif
