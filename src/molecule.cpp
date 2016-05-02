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

      
      double val = (double) dis(gen);
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


double Molecule::determinant(std::vector<std::vector<double>> A, int N){

  double det = 1.0;

  for(int rc=0; rc<N; rc++){
    if (abs(A[rc][rc]) < err){
      int new_col = -1;
      for(int i=rc+1; i < N; i++){
        if (abs(A[i][rc]) > err){
          new_col = i;
          break;
        }
      }

      if (new_col != -1){
        det = det * (-1.0);
        for(int i=0;i<N; i++){
          double temp = A[rc][i];
          A[rc][i] = A[new_col][i];
          A[new_col][i] = temp;
        }
      }

    }
    for(int row_below=rc+1; row_below<N; row_below++){
      if (abs(A[row_below][rc]) > err){
        double val = fmod(A[row_below][rc]/A[rc][rc], prime);
        for(int i=0; i<N; i++)
          A[row_below][i] = fmod(A[row_below][i] - val*A[rc][i], prime);
      }
    }

  }

  for(int i=0; i<N; i++)
    det = fmod(det*A[i][i], prime);

  return det;
}


/* inverse - Computes the inverse of an NxN matrix  A in place
*/
void Molecule::inverse(std::vector<std::vector<double>> &C, int N, std::map<int, int> excl){
  

  double err_power = err*(((double) 4*excl.size())/((double) N)+10.0);
  err_power = .0000000000001;
  printf("err_power: %.3e\n", err_power);



  int j = 0;
  while (j<N){
    
    if (excl.count(j)>0){
      j += 1;
      continue;
    }


    if (j % 10 == 0){
      for (int r=j+1; r<N;r++){
        for(int c=0; c<N; c++){
          C[r][c]  = C[r][c]*1000;
        }
      }
    }


    if (std::abs(C[j][j]) < err){
      int k = -1;

      for(int new_col=j+1; new_col<N; new_col++){
        if (excl.count(new_col) == 0 && std::abs(C[new_col][j]) > 0){
          k = new_col;
          break;
        }
      }

      if (k==-1){
        printf("k is -1 \n");
        printf("j: %d\n", j);
        for(int new_col=j+1; new_col<N; new_col++){
          printf("in excl: %s value: %.4e above err: %s\n", excl.count(new_col) == 0 ? "true" : "false", 
            C[new_col][j],
            std::abs(C[new_col][j]) > err ? "true" : "false");
        }
      }

      for(int row=0; row<2*N; row++){
        if (excl.count(row) == 0){
          C[j][row] = C[j][row] + C[k][row],prime;
        }
      }
    }

    double ajj = C[j][j];

    for(int row=0; row<2*N; row++){
      if (excl.count(row) == 0)
        C[j][row] = C[j][row]/ajj, prime;
    }

    for(int col=0; col < N; col++){
      if ((col !=j) && excl.count(col)==0){

        double aij = C[col][j];

        for(int row=0;row<2*N;row++){
          if (excl.count(row) ==0){
            C[col][row] = C[col][row] - aij*C[j][row], prime;
          }
        }
      }
    }
    j++;
  }
}

std::vector<std::vector<double>> Molecule::copy_graph(){

  int N = numberOfAtoms();
  vector<vector<double>> copy = vector<vector<double>>(N);
  for(int i=0; i<N; i++)
    copy[i] =  vector<double>(N);

  for (int i=0; i< N; i++)
    for (int j=0; j<N; j++)
        copy[i][j] = graph[i][j];
  
  return copy;
}

std::vector<std::vector<double>> Molecule::copy_augC(){

  int N = numberOfAtoms();
  vector<vector<double>> copy = vector<vector<double>>(N);
  for(int i=0; i<N; i++)
    copy[i] = vector<double>(2*N);

  for (int i=0; i< N; i++)
    for (int j=0; j<2*N; j++)
        copy[i][j] = augC[i][j];
  
  return copy;
}


std::vector<std::tuple<int, int>> Molecule::matching(){

  unsigned int n = numberOfAtoms();
  std::vector<std::vector<double>> graph_copy = copy_graph();
  double det = determinant(graph_copy, n);

  printf("determinant: %f\n", det);
  printf("det condition: %s\n", (std::abs(det) < err) ? "true" : "false");

  if (std::abs(det) < err){
    return std::vector<std::tuple<int, int>>();
  }


  printf("current error is: %.3e\n", err);

  std::vector<std::tuple<int, int>> M = std::vector<std::tuple<int, int>>();
  std::map<int, int> excl = std::map<int, int>();

  while (M.size() < n/2){
    
    std::vector<std::vector<double>> N = copy_augC();
    inverse(N, n, excl);
    int first_row = -1;

    for(int row=0; row<n; row++){
      if (excl.count(row) == 0){
        first_row = row;
        break;
      }
    }

    if (first_row == -1){
      printf("FIRST ROW is -1");
      printf("matching size = %d\n", M.size());
      return std::vector<std::tuple<int, int>>();
    }
  

    int j_col = -1;
    for(int col=n; col<2*n; col++){
      if (excl.count(col-n) == 0 && abs(N[first_row][col])>err && graph[first_row][col-n] != 0){
        j_col = col - n;
        break;
      }
    }

    if (j_col == -1){
      printf("j_col is -1\n");
      printf("matching size = %d\n", M.size());
      return std::vector<std::tuple<int, int>>();
    }


    printf("size of matching: %d\n", M.size());
    M.push_back(std::make_tuple(first_row, j_col));

    //printf("(%d, %d)\n", first_row, j_col);
    //printf("Excl values: ");
    //for(auto it= excl.begin(); it != excl.end(); ++it){
      //printf(" %d ", it->first);
    //}
    //printf("\n");
    //printMatrix(N);   
    excl.insert(std::pair<int,int>(first_row, 1));
    excl.insert(std::pair<int,int>(j_col, 1));
  }
  return M;
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

void Molecule::printMatrix(std::vector<std::vector<double>> A){
  int count = 0;

  for(auto it= A.begin(); it != A.end(); ++it){
    printf(" %d ---->  ", count);
    for(auto f= it->begin(); f != it->end(); ++f){
      printf(" %f ", *f);

    }
    printf("\n");
    count++;
  }
}



void Molecule::printGraph(){
  int count = 0;

  for(auto it= graph.begin(); it != graph.end(); ++it){
    printf("[");
    for(auto f= it->begin(); f != it->end(); ++f){
      printf("%f,", *f);

    }
    printf("],\n");
    count++;
  }
  printf("%d\n",numBonds);

}


void Molecule::printDeterminant(){

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
  graph = vector<vector<double>>(N);
  augC = vector<vector<double>>(N);

  for(int i=0; i<N; ++i){
    graph[i] =  vector<double>(N);
    augC[i] = vector<double>(2*N);
  }

  for (int j=0; j< N; j++){
    for (int k=0; k<N; k++){
      if (j==k){
        augC[j][k+N] = 1;
      }
    }

  }

  
}

