// Simple Molecule class

#include "molecule.h"
#include <random>
#include <iostream>
//#include <cuda_runtime.h>
//#include <cuda.h>


using namespace std;
using namespace Eigen;

extern std::vector<std::tuple<int, int>> setup(double *cudaGraph, vector<vector<double>> host_graph, int N, double  err);

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

  //#pragma parallel for schedule(static)
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

      cudaGraph[atom->id()*n + nbr->id()] = val;
      cudaGraph[nbr->id()*n + atom->id()] = -1.0*val;
      numedges++;

      addBond(atom,nbr,1);
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

    #pragma omp parallel for schedule(static)
    for (int col=k; col<m; col++){
        double temp = U(k,col);
        U(k,col) = U(i,col);
        U(i,col) = temp;
    }

    #pragma omp parallel for schedule(static)
    for (int col=0; col<k; col++){
        double temp = L(k,col);
        L(k,col) = L(i,col);
        L(i,col) = temp;
    }

    #pragma omp parallel for schedule(static)
    for (int col=0; col<m; col++){
        double temp = P(k,col);
        P(k,col) = P(i,col);
        P(i,col) = temp;
    }

    #pragma omp parallel for schedule(static)
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
    #pragma omp parallel for reduction(-:rhs)
    for (int j=0; j<i; j++){
      rhs -= L(i,j)*y(j);
    }
    y(i) = rhs;
  }

  for (int i=m-1; i>=0; i--){
    double rhs = y(i);

    #pragma omp parallel for reduction(-:rhs)
    for (int j=m-1; j> i; j--){
      rhs -= U(i,j)*x(j);
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


std::vector<std::tuple<int, int>> Molecule::CUDAMatching(){
  unsigned int n = numberOfAtoms();
  printf("Error in mol.cpp: %.3e\n", err);
  std::vector<std::tuple<int, int>> M = setup(cudaGraph, graph, n, err);
  printMatching(M);
  return M;
}

std::vector<std::tuple<int, int>> Molecule::matching(){

  unsigned int n = numberOfAtoms();
  int matrix_size = n;
  std::vector<std::tuple<int, int>> M = std::vector<std::tuple<int, int>>();


  MatrixXd A(matrix_size, matrix_size);
  VectorXd b = VectorXd::Zero(matrix_size);
  b[0] = 1;

  //Set A
 #pragma omp parallel for schedule(static)
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
      VectorXd x = solve(A,b,matrix_size);
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
          int counter_col=0;
          for (int j=0; j<matrix_size; j++){
            if (j!=0 && j!=row_j){
              //A(row_counter, col_count[j]-1) = A_copy(i,j);
              A(row_counter, counter_col) = A_copy(i,j);
              counter_col++;
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
  return M;
}


void Molecule::addMatchedBonds(std::vector<std::tuple<int, int>> M){

    int sizeMatching = M.size();
    std::map<int, int> matchingMap;
    for(int i=0; i<sizeMatching; i++){
        int a = std::get<0>(M[i]);
        int b = std::get<1>(M[i]);
        matchingMap.insert(std::pair<int,int>(a,b));
    }


    for(int j=0; j<_bonds.size(); j++){

        int id_first = getfirst(j);
        int id_second= getsecond(j);

        int min = std::min(id_first,id_second);
        int max = std::max(id_first, id_second);

        if (matchingMap.count(min)==1 &&
            matchingMap[min]==max){
            update(j);
        }
    }

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
//  cudaMallocHost((void**)cudaGraph,N*N);

  cudaGraph = (double *)calloc(N*N, sizeof(double));

  //#pragma parallel for schedule(static)
  for(int i=0; i<N; ++i){
    graph[i] =  vector<double>(N);
  }


}

