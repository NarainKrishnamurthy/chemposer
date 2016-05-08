#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include "molecule.h"
#include <vector>

using namespace std;

#define THREADS_PB 256

__device__ double device_graph;
__device__ double *A;
__device__ double *P;
__device__ double *L;
__device__ double *U;
__device__ double *Pb;
__device__ double *b;
__device__ double *x;
__device__ double *y;
__device__ int *M1;
__device__ int *M2;


__global__ void kernelRowSwap(double *M, int start, int end, int k, int i, int N){
    int threadIndex = blockDim.x*blockIdx.x + threadIdx.x;

    if (threadIndex >= (start - end))
        return;

    double temp = M[k*N + start + threadIndex];
    M[k*N + start + threadIndex] = M[i*N + start + threadIndex];
    M[i*N + start + threadIndex] = temp;
}



__global__ void kernelSetLURow(int k, int N){
    int j = blockDim.x*blockIdx.x + threadIdx.x;

    if (j >= N || j<k+1)
        return;

    L[j*N + k] = U [j*N + k]/U[k*N+k];
    for (int col=k; col<N; col++){
        U[j*N + col] -= L[j*N+k]*U[k*N+col];
    }
}


__global__ void kernelMakeIdentity(double *I, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    I[i*N+j] = (i==j) ? 1.0 : 0.0;
}

__global__ void copyAtoU(int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    U[i*N+j] = A[i*N+j];
}


//REQUIRES: row_counter < i
__global__ void kernelresizeA(int row_counter, int i, int row_j, int N){
    int new_size = N-2;
    int j = blockDim.x*blockIdx.x + threadIdx.x;

    if (j >= N || j == 0 || j == row_j)
        return;

    if (j < row_j){
        A[row_counter*new_size + j - 1] = A[i*N + j];
    }
    else{
        A[row_counter*new_size + j - 2] = A[i*N + j];
    }
}


__global__ void kernelSetPb(int N){
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i >= N)
        return;

    if (P[i*N] == 1.0){
        Pb[i] = 1.0;
    } else{
        Pb[i] = 0.0;
    }
}

__global__ void kernelSetb(int N){
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i >= N)
        return;

    b[i] = i == 0 ? 1.0 : 0.0;
}


__global__ void  kernelSequentialSolve(int N){


  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      U[i*N + j] = A[i*N + j];
  
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (i==j){
        L[i*N + j] = 1.0;
        P[i*N + j] = 1.0;
      } else {
        L[i*N + j] = 0.0;
        P[i*N + j] = 0.0;
      }
    }
  }

  for (int k=0; k<N-1; k++){
    int i = -1;
    double max = -1;

    for (int row = k; row<N; row++){
      if (abs(U[row*N + k]) > max){
        i = row;
        max = abs(U[row*N+k]);
      }
    }

    for (int col=k; col<N; col++){
        double temp = U[k*N+col];
        U[k*N+col] = U[i*N+col];
        U[i*N+col] = temp;
    }

    for (int col=0; col<k; col++){
        double temp = L[k*N+col];
        L[k*N+col] = L[i*N+col];
        L[i*N+col] = temp;
    }

    for (int col=0; col<N; col++){
        double temp = P[k*N+col];
        P[k*N+col] = P[i*N+col];
        P[i*N+col] = temp;
    }

    for (int j=k+1; j<N; j++){
      L[j*N+k] = U[j*N+k]/U[k*N+k];
      for (int col=k; col<N; col++){
        U[j*N+col] = U[j*N+col] - L[j*N+k]*U[k*N+col];
      }
    }
  }

  for (int i=0; i<N; i++){
    x[i] = 0.0;
    y[i] = 0.0;
  }

  for (int i=0; i<N; i++){
    double sum = 0.0;
    for (int j=0; j<N;j++){
      sum += P[i*N+j]*b[j];
    }
    Pb[i] = sum;
  }

  for (int i = 0; i < N; i++){
    double rhs = Pb[i];

    for (int j=0; j<i; j++){
      rhs -= L[i*N+j]*y[j];
    }
    y[i] = rhs;
  }

  for (int i=N-1; i>=0; i--){
    double rhs = y[i];

    for (int j=N-1; j> i; j--){
      rhs -= U[i*N+j]*x[j];
    }
    x[i] = rhs/U[i*N+i];
  }
}



std::vector<std::tuple<int, int>> matching(std::vector<std::vector<double>> graph, int n, double err){

  int matrix_size = n;
  std::vector<std::tuple<int, int>> M = std::vector<std::tuple<int, int>>();

  vector<int> rc_map = vector<int>();
  for (int i=0; i<n; i++){
    rc_map.push_back(i);
  }

  
  double* host_x = (double *)calloc(n, sizeof(double));

  while (M.size() < n/2){
      dim3 seqblockDim(1024,1,1);
      dim3 seqgridDim((1 + seqblockDim.x - 1)/seqblockDim.x,1,1);

      kernelSequentialSolve<<<seqgridDim,seqblockDim>>>(matrix_size);
      cudaThreadSynchronize();
      int row_j = -1;
      
      //COpy x to host

      cudaMemcpy(host_x, x, matrix_size*sizeof(double), cudaMemcpyDeviceToHost);

      printf("returned x\n");
      for (int i=0; i<matrix_size; i++)
        printf("%f\n", host_x[i]);

      for (int row=0; row< matrix_size; row++){
        int true_first_col = rc_map[0];
        int true_row = rc_map[row];
        if (abs(host_x[row]) > err && graph[true_row][true_first_col]!=0){
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
      int row_counter = 0;


      for (int i=0; i<matrix_size; i++){
        if (i != 0 && i != row_j){

          dim3 blockDim(1024,1,1);
          dim3 gridDim((matrix_size + blockDim.x - 1)/blockDim.x,1,1);

          kernelresizeA<<<gridDim, blockDim>>>(row_counter, i, row_j, matrix_size);
          cudaThreadSynchronize();
          row_counter++;
        }
      }

      matrix_size -= 2;
      rc_map.erase(rc_map.begin() + row_j);
      rc_map.erase(rc_map.begin());
  }
  return M;
}


std::vector<std::tuple<int, int>> setup(std::vector<std::vector<double>> host_graph, int N, int err){
    
    dim3 blockDim(1024,1,1);
    dim3 gridDim((N + blockDim.x - 1)/blockDim.x,1,1);


    cudaMalloc((void**)&A, sizeof(double)*N*N);
    cudaMalloc((void**)&P, sizeof(double)*N*N);
    cudaMalloc((void**)&L, sizeof(double)*N*N);
    cudaMalloc((void**)&U, sizeof(double)*N*N);
    cudaMalloc((void**)&Pb, sizeof(double)*N);

    for(int i=0; i<N; i++){
      cudaMemcpy(A + i*N*sizeof(double), &host_graph[i], N*sizeof(double), cudaMemcpyHostToDevice); 
    }

    kernelSetb<<<gridDim, blockDim>>>(N);
    cudaThreadSynchronize();

    return matching(host_graph, N, err);  
}

/*
    kernelRowSwap<<<gridDim, blockDim>>>(A,start,end,k,i,N);
    kernelSetLURow<<<gridDim, blockDim>>>(k,  N);    
    kernelMakeIdentity<<<gridDim, blockDim>>>(L, N);
    copyAtoU<<<gridDim,blockDim>>>(N);

    int row_j = 0;
    int row_counter = 0;
    kernelresizeA<<<gridDim, blockDim>>>(row_counter, i, row_j,  N);
    kernelSetPb<<<gridDim, blockDim>>>(N);
*/

/*
__global__ void kernelForwardSolve(double *y, double *L, double* Pb, int N){
    for (int i = 0; i < N; i++){
        double rhs = Pb[i];

        for (int j=0; j<i; j++){
          rhs -= L[i*N +j]*y[j];
        }
        y[i] = rhs;
    }
}


__global__ void kernelBackwardSolve(double *x, double *U, double* y int N){
      for (int i=N-1; i>=0; i--){
        double rhs = y[i];

        for (int j=N-1; j> i; j--){
          rhs -= U[i*N+j]*x[j];
        }
        x[i] = rhs/U[i*N+i];
      }
}

void solve(int N){


    //A, U, L, P, y, x, Pb

    for (int k=0; k<N-1; k++){
        int i = -1;
        double max = -1.0;

        for (int row = k; row<N; row++){
          if (abs(U(row, k)) > max){
            i = row;
            max = abs(U(row, k));
          }
        }

        kernelRowSwap(U, k, N, k,i, N);
        kernelRowSwap(L, 0, k, k,i, N);
        kernelRowSwap(P, 0, N, k,i, N);

        kernelSetLURow(L,U,k,N);
    }

    kernelSetPb(Pb, P, N);
    kernelForwardSolve(y,L,Pb,N);
    kernelBackwardSolve(x,U,y,N);

    //x will now be set to the valid x
}

*/



