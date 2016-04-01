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
 { Atom *a = new Atom(element, x, y, z); _atoms.push_back(a); }
 void addBond(Bond *bond) { _bonds.push_back(bond); }

 unsigned int numberOfAtoms() { return _atoms.size(); }
 unsigned int numberOfBond() { return _bonds.size(); }

 void clear() { _atoms.clear(); _bonds.clear(); }

 protected:
  std::vector<Atom*> _atoms;
  std::vector<Bond*> _bonds;
};

#endif
