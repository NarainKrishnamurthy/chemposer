// Simple Molecule class

#include "molecule.h"
#include <random>

using namespace std;

bool sortAtomZ(const pair<Atom*,double> &a, const pair<Atom*,double> &b)
{   return (a.second < b.second); }

void Molecule::perceiveBonds() {
  // find the initial bonded atoms
  unsigned int max;
  double maxrad = 0.0;
  bool unset = false;
  Atom *atom,*nbr;
  vector<Atom*>::iterator i;
  vector<pair<Atom*,double> > zsortedAtoms;
  vector<double> rad;

  rad.resize(_atoms.size());

  for (i = _atoms.begin(); i < _atoms.end(); ++i) {
    pair<Atom*,double> entry(*i, (*i)->z());
    zsortedAtoms.push_back(entry);
  }
  sort(zsortedAtoms.begin(), zsortedAtoms.end(), sortAtomZ);
  max = zsortedAtoms.size();

  // find the maximum radius for the atoms
  for (unsigned int j = 0; j < max; ++j) {
    atom   = zsortedAtoms[j].first;
    rad[j] = atom->radius();
    maxrad = std::max(rad[j],maxrad);
  }

  int numedges = 0;
  int idx1, idx2;
  double d2,cutoff,zd;

  int n = numberOfAtoms();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(1, n*n);

  for (unsigned int j = 0; j < max ; ++j) {
    double maxcutoff = SQUARE(rad[j]+maxrad+0.45);
    atom = zsortedAtoms[j].first;

    for (unsigned int k = j + 1 ; k < max ; k++ ) {
      nbr = zsortedAtoms[k].first;

      // bonded if closer than elemental Rcov + tolerance
      cutoff = SQUARE(rad[j] + rad[k] + 0.45);
      // current z distance
      zd  = SQUARE(atom->z() - nbr->z());
      // bigger than max cutoff, which is determined using largest radius,
      // not the radius of k (which might be small, ie H, and cause an early  termination)
      // since we sort by z, anything beyond k will also fail
      if (zd > maxcutoff )
        break;

      d2  = SQUARE(atom->x() - nbr->x());
      if (d2 > cutoff)
        continue; // x's bigger than cutoff
      d2 += SQUARE(atom->y() - nbr->y());
      if (d2 > cutoff)
        continue; // x^2 + y^2 bigger than cutoff
      d2 += zd;
      if (d2 > cutoff)
        continue; // too far
      if (d2 < 0.16) // 0.4 * 0.4 = 0.16
        continue; // too close

      
      float val = (float) dis(gen);
      if (atom->id() > nbr->id()){
         val = -1.0 * val;
      } 
      graph[atom->id()][nbr->id()] = val;
      graph[nbr->id()][atom->id()] = -1.0*val;

      augC[atom->id()][nbr->id()] = val;
      augC[nbr->id()][atom->id()] = -1.0*val;

      numedges++;
    } // end inner loop
  } // end outer loop
  setNumBonds(numedges);
  //printGraph();

} // end perceiveBonds

/*

int determinant(N, A){
    int det = 1;
    for (int i=0; i<N; i++){
        int aii = A[i][i];
        det = det*aii;
        //#pragma omp parallel for
        for (int j=i+1; j<N; j++){
            int z = A[j][i]/aii;
            for (int k=i+1; k<N; k++){
                A[j][k] = A[j][k] - z*A[i][k];
            }

        }
    }
    return det;
}

*/

/* inverse - Computes the inverse of an NxN matrix  A in place
*/
void  Molecule::inverse(std::map<int, int> *excl_rows, std::map<int, int> *excl_cols){
  int N = numberOfAtoms();
  int j = 0;

  while (j<N){

    if (excl_cols->count(j)>0){
      j += 1;
      continue;
    }
    if (augC[j][j] == 0){
      int k = -1;

      for(int new_col=j+1; new_col<N; new_col++){
        if (excl_cols->count(new_col) == 0 && augC[new_col][j] != 0){
          k = new_col;
          break;
        }
      }

      if (k==-1){
        printf("k is -1 \n");
      }

      for(int row=0; row<2*N; row++){
        if (excl_rows->count(row) == 0){
          augC[j][row] += augC[k][row];
        }
      }

      
    }

    int ajj = augC[j][j];

    for(int row=0; row<2*N; row++){
      if (excl_rows->count(row) == 0)
        augC[j][row] /= ajj;
    }

    for(int col=0; col < N; col++){
      if ((col !=j) && excl_cols->count(col)==0){

        int aij = augC[col][j];

        for(int row=0;row<2*N;row++){
          if (excl_rows->count(row) ==0){
            augC[col][row] = augC[col][row] - aij*augC[j][row];
          }
        }
      }


    }

    j++;

  }

}


void Molecule::doMatching() {

  // loop through the bonds
  vector<Bond*>::iterator it;
  for (it = _bonds.begin(); it < _bonds.end(); ++it) {
    // do something here
  } // end bond loop

}


void Molecule::printMolecule(){

 for(auto it= _atoms.begin(); it != _atoms.end(); ++it){
    printf(" %d \n", (*it)->id());

  }

}

void Molecule::printGraph(){
  int count = 0;

  for(auto it= graph.begin(); it != graph.end(); ++it){
    printf(" %d ---->  ", count);
    for(auto f= it->begin(); f != it->end(); ++f){
      printf(" %f ", *f);

    }
    printf("\n");
    count++;
  }
  printf("%d\n",numBonds);

}


void Molecule::printAugMatrix(){
  int count = 0;

  for(auto it= augC.begin(); it != augC.end(); ++it){
    printf(" %d ---->  ", count);
    for(auto f= it->begin(); f != it->end(); ++f){
      printf(" %.1lf ", *f);

    }
    printf("\n");
    count++;
  }
 

}

void Molecule::initializeGraph(){

  int N = numberOfAtoms();
  graph = vector<vector<float>>(N);
  augC = vector<vector<float>>(N);

  for(int i=0; i<N; ++i){
    graph[i] =  vector<float>(N);
    augC[i] = vector<float>(2*N);
  }

  for (int j=0; j< N; j++){
    for (int k=0; k<N; k++){
      if (j==k){
        augC[j][k+N] = 1;
      }
    }

  }

  
}

