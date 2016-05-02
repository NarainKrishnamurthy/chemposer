// Simple Molecule class

#include "molecule.h"
#include <random>
#include <iostream>

using namespace std;
using namespace Eigen;

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


VectorXd  Molecule::solve(MatrixXd A,  VectorXd b, int m){
  MatrixXd U = MatrixXd(A);
  MatrixXd L = MatrixXd::Identity(m, m);
  MatrixXd P = MatrixXd::Identity(m, m);
  for (int k=0; k<m-1; k++){
    int i = -1;
    double max = -1;
    for (int row = k; row<m; row++){
      if (abs(U(row, k)) > max){
        i = row;
        max = abs(U(row, k));
      }
    }

    for (int col=k; col<m; col++){
        double temp = U(k,col);
        U(k,col) = U(i,col);
        U(i,col) = temp;
    }

    for (int col=0; col<k; col++){
        double temp = L(k,col);
        L(k,col) = L(i,col);
        L(i,col) = temp;
    }

    for (int col=0; col<m; col++){
        double temp = P(k,col);
        P(k,col) = P(i,col);
        P(i,col) = temp;
    }

    for (int j=k+1; j<m; j++){
      L(j,k) = U(j,k)/U(k,k);
      for (int col=k; col<m; col++){
        U(j,col) = U(j,col) - L(j,k)*U(k,col);
      }
    }
  }

  VectorXd x = VectorXd::Zero(m);
  VectorXd y = VectorXd::Zero(m);
  VectorXd Pb = P*b;



  for (int i = 0; i < m; i++){
    double rhs = Pb(i);
    for (int j=0; j<i; j++){
      rhs = rhs - L(i,j)*y(j);
    }
    y(i) = rhs;
  }

  for (int i=m-1; i>=0; i--){
    double rhs = y(i);
    for (int j=m-1; j> i; j--){
      rhs = rhs - U(i,j)*x(j);
    }
    x(i) = rhs/U(i,i);
  }

  //const IOFormat fmt(2, DontAlignCols, "\t", " ", "", "", "", "");
  //cout << (U*x).format(fmt) << endl;
  //cout << (y).format(fmt) << "\n"<< endl;

  //bool is_approx = (U*x).isApprox(y, 1.0);
  //printf("Matrices are %s\n", is_approx ? "SIMILAR" : "NOT SIMILAR");


  //const IOFormat fmt(2, DontAlignCols, "\t", " ", "", "", "", "");
  //cout << ((P*A) - (L*U)).format(fmt) << "\n"<< endl;
  //bool is_approx = (P*A).isApprox(L*U, 1.0);
  //printf("Matrices are %s\n", is_approx ? "SIMILAR" : "NOT SIMILAR");
  return x;
}


std::vector<std::tuple<int, int>> Molecule::matching(){

  unsigned int n = numberOfAtoms();
  int matrix_size = n;
  std::vector<std::tuple<int, int>> M = std::vector<std::tuple<int, int>>();


  MatrixXd A(matrix_size, matrix_size);
  VectorXd b = VectorXd::Zero(matrix_size);
  b[0] = 1;

  //Set A
  for (int row=0; row<n; row++){
    for (int col=0; col<n; col++){
      A(row, col) = .1*graph[row][col];
    }
  }

  vector<int> rc_map = vector<int>();
  for (int i=0; i<n; i++){
    rc_map.push_back(i);
  }

  while (M.size() < n/2){
      solve(A,b,matrix_size);
      VectorXd x = solve(A,b,matrix_size);
      //VectorXd x = A.partialPivLu().solve(b);

      int row_j = -1;
      
      for (int row=0; row< matrix_size; row++){
        int true_first_col = rc_map[0];
        int true_row = rc_map[row];
        if (abs(x[row]) > err && graph[true_row][true_first_col]!=0){
          row_j = row;
          break;
        }
      }

      if (row_j == -1){
        printf("row_j is -1\n");
        printf("matching size = %d\n", M.size());
        return std::vector<std::tuple<int, int>>();
      }

      M.push_back(std::make_tuple(rc_map[0], rc_map[row_j]));
      MatrixXd A_copy = MatrixXd(A);

      int row_counter = 0;
      for (int i=0; i<matrix_size; i++){
        if (i != 0 && i != row_j){
          int col_counter = 0;
          for (int j=0; j<matrix_size; j++){
            if (j!=0 && j!=row_j){
              A(row_counter, col_counter) = A_copy(i,j);
              col_counter++;
            }
          }
          row_counter++;
        }
      }
      matrix_size -= 2;
      A.conservativeResize(matrix_size, matrix_size);
      b.conservativeResize(matrix_size,1);
      if (matrix_size != 0){
        b[0] = 1;
      }

      rc_map.erase(rc_map.begin() + row_j);
      rc_map.erase(rc_map.begin());

  }
  printMatching(M);
  return M;
}

void Molecule::printMatching(std::vector<std::tuple<int, int>> M){
  for(auto it= M.begin(); it != M.end(); ++it){
    printf("(%d, %d)\n", get<0>(*it), get<1>(*it));
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

/*
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


*/
