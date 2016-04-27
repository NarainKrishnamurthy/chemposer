// Simple Molecule class

#include "molecule.h"

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

      unsigned int n = numberOfAtoms();
      graph[atom->id()][nbr->id()] = 1;
      graph[nbr->id()][atom->id()] = 1;

      augC[atom->id()][n+nbr->id()] = 1;
      augC[nbr->id()][n+atom->id()] = 1;

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
void  Molecule::inverse(){
    int N = numberOfAtoms();
    int j = 0;
    while (j < N) {
        int k = 0;
        //IF... add a non-zero row to row j
        if (augC[j][j] == 0){
            //Find the non-zero k
            //#pragma omp parallel for
            for (int i=0; i<N; i++){
                if (augC[k][j]!=0){
                    k = i;
                    break;
                }
            }
            //Add the row k to row j
            //#pragma omp parallel for
            for (int i=j; i<N+j; i++){
                augC[j][i] = augC[j][i] + augC[k][i];
            }
        }
        float ajj = augC[j][j];
        
        //Divide out ajj
        //#pragma omp parallel for
        for (int i=j; i<N+j; i++){
            augC[j][i] = augC[j][i]/ajj;
        }

        //Subtract row j multiplied by appropriate constant from other rows
        for (int i=0; i<N;i++){
            if (i!=j){
                float aij = augC[i][j];
                for (int r=j; r<N+j; r++){
                    augC[i][r] = augC[i][r] - aij*augC[j][r-j];
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

  graph = vector<vector<float>>(numberOfAtoms());
  augC = vector<vector<float>>(numberOfAtoms());

  for(int i=0; i<numberOfAtoms(); ++i){
    graph[i] =  vector<float>(numberOfAtoms());
    augC[i] = vector<float>(2*numberOfAtoms());
  }

  for (int j=0; j< numberOfAtoms(); j++){
    for (int k=0; k<numberOfAtoms(); k++){
      if (j==k){
        augC[j][k] = 1;
      }

    }

  }

  
}

