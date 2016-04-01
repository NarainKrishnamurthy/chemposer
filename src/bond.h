// Simple Bond class

#ifndef BOND_H
#define BOND_H

class Atom;

class Bond {
 public:
 Bond(): _start(NULL), _end(NULL) {};
 Bond(Atom *a, Atom *b): _start(a), _end(b) {};
 ~Bond() {}

 void setStart(Atom *a) { _start = a; }
 Atom *start() { return _start; }

 void setEnd(Atom *b) { _start = b; }
 Atom *end() { return _end; }

 void set(Atom *a, Atom *b) { _start = a; _end = b; }

 Atom *neighbor(Atom *a)
 { if (a == _start) return _end;
   else if (a == _end) return _start;
   else return NULL;
 }

 protected:
  Atom* _start;
  Atom* _end;
};

#endif
