// Simple Bond class

#ifndef BOND_H
#define BOND_H

#include <cstddef>
class Atom;

class Bond {
 public:
 Bond(): _start(NULL), _end(NULL), _type(0) {};
 Bond(Atom *a, Atom *b): _start(a), _end(b), _type(0) {};
 ~Bond() {}

 void setStart(Atom *a) { _start = a; }
 Atom *start() { return _start; }

 void setEnd(Atom *b) { _start = b; }
 Atom *end() { return _end; }

 void set(Atom *a, Atom *b) { _start = a; _end = b; }

 Atom *neighbor(Atom *a)
 { if (a == _start) return _end;
   else if (a == _end) return _start;
   return NULL;
 }

 int firstID(){ return _start->id();}

 int secondID(){ return _end->id();}

 int getType() {return _type;}

 unsigned short type() { return _type; }
 void setType(unsigned short t) { _type = t; }

 protected:
  Atom* _start;
  Atom* _end;

  unsigned short _type; // 0 = unknown, 1, 2, 3, etc.
};

#endif
