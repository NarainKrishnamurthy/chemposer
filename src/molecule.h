// Simple Molecule class

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>

#include "atom.h"
#include "bond.h"

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
 unsigned int numberOfBonds() { return _bonds.size(); }

 std::vector<Atom*> atoms() { return _atoms; }
 std::vector<Bond*> bonds() { return _bonds; }

 void clear() { _atoms.clear(); _bonds.clear(); }

 void perceiveBonds();
 void doMatching();

 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;
};

#endif
