#include "atom.h"
#include "bond.h"

using namespace std;

bool Atom::isBonded(Atom &a)
{
  bool bonded = false;

  // loop through the bonds
  vector<Bond*>::iterator it;
  for (it = _bonds.begin(); it < _bonds.end(); ++it) {
    if ((*it)->neighbor(this) == &a) {
      bonded = true;
      break;
    }
  }

  return bonded;
}
