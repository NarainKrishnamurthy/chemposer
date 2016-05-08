#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include "molecule.h"
#include <vector>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

using namespace std;

#define THREADS_PB 256

__device__ double *A;
__device__ double *L;
__device__ double *U;
__device__ double *P;
__device__ double *Pb;
__device__ double *b;
__device__ double *x;
__device__ double *y;


__global__ void kernelRowSwap(double *M, int start, int end, int k, int* i_ptr, int N){
    
    int i = *i_ptr;
    int col = blockDim.x*blockIdx.x + threadIdx.x;

    if (col < start || col >= end)
        return;

    double temp = M[k*N + col];
    M[k*N + col] = M[i*N + col];
    M[i*N + col] = temp;
}



__global__ void kernelSetLU(double *L, double *U, int k, int N){
    int j = blockDim.x*blockIdx.x + threadIdx.x;

    if (j >= N || j<k+1)
        return;

    L[j*N + k] = U [j*N + k]/U[k*N+k];
    for (int col=k; col<N; col++){
        U[j*N + col] -= L[j*N+k]*U[k*N+col];
    }
}


__global__ void kernelInitLP(double *L, double *P, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    L[i*N+j] = (i==j) ? 1.0 : 0.0;
    P[i*N+j] = (i==j) ? 1.0 : 0.0;
}

__global__ void kernelcopyAtoU(double *A, double *U, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    U[i*N+j] = A[i*N+j];
}


//REQUIRES: row_counter < i
__global__ void kernelresizeA(double *A, int row_counter, int i, int row_j, int N){
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


__global__ void kernelSetPb(double *Pb, double *P, int N){
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i >= N)
        return;

    if (P[i*N] == 1.0){
        Pb[i] = 1.0;
    } else{
        Pb[i] = 0.0;
    }
}

__global__ void kernelSetb(double *b, int N){
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i >= N)
        return;

    b[i] = (i==0) ? 1.0 : 0.0;
}

__global__ void kernelResetxy(double *x, double *y, int N){
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i >= N)
        return;

    x[i] = 0.0;
    y[i] = 0.0;
}

__global__ void kernelCopyRowA(double *A, double* row, int i, int N){
  int j = blockDim.x*blockIdx.x + threadIdx.x;
  if (j>= N)
    return;

  A[i*N + j] = row[j];
}

__global__ void kernelScaleA(double *A, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    A[i*N+j] = A[i*N+j];
}
__global__ void kernelGetMaxUk(double *U, int k, int* i, int N){
  double max = -1.0;
  *i = -1;
  for (int row = k; row<N; row++){
    if (abs(U[row*N+k]) > max){
      *i = row;
      max = abs(U[row*N+k]);
    }
  }
}

__global__ void kernelForwardMap(double *temp, double *L, double *y, int i, int N){
  int j = blockDim.x*blockIdx.x + threadIdx.x;
  if (j>= i)
    return;

  temp[j] = L[i*N+j]*y[j];
}

__global__ void kernelForwardSolve(double* Pb, double *y, double result, int i, int N){

  double rhs = Pb[i];
  rhs -= result;
  y[i] = rhs;
}


__global__ void kernelSequentialHelpAfter(double *A, double *L, double *U, double *P,
  double *Pb, double *b, double *x, double *y, int N){

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

void solve(int N){
  
  dim3 blockDimVector(1024,1,1);
  dim3 gridDimVector((N + blockDimVector.x - 1)/blockDimVector.x,1,1);

  dim3 blockDimArray(1024,1,1);
  dim3 gridDimArray((N + blockDimArray.x - 1)/blockDimArray.x,1,1);

  blockDimArray.x = 32;
  blockDimArray.y = 32;
  gridDimArray.x = (N + blockDimArray.x -1)/ blockDimArray.x;
  gridDimArray.y = (N + blockDimArray.y -1)/ blockDimArray.y;

  kernelcopyAtoU<<<gridDimArray, blockDimArray>>>(A,U,N);
  cudaThreadSynchronize();

  kernelInitLP<<<gridDimArray, blockDimArray>>>(L,P,N);
  cudaThreadSynchronize();

  for (int k=0; k<N-1; k++){
    int *i;
    cudaMalloc((void**)&i, sizeof(int));

    kernelGetMaxUk<<<1, 1>>>(U,k,i,N);
    cudaThreadSynchronize();
    
    kernelRowSwap<<<gridDimVector, blockDimVector>>>(U, k, N, k,i, N);
    cudaThreadSynchronize();

    kernelRowSwap<<<gridDimVector, blockDimVector>>>(L, 0, k, k,i, N);
    cudaThreadSynchronize();

    kernelRowSwap<<<gridDimVector, blockDimVector>>>(P, 0, N, k,i, N);
    cudaThreadSynchronize();

    cudaFree(i);
    kernelSetLU<<<gridDimVector, blockDimVector>>>(L,U,k,N);
    cudaThreadSynchronize();
  }

  kernelResetxy<<<gridDimVector, blockDimVector>>>(x,y,N);
  cudaThreadSynchronize();

  kernelSetPb<<<gridDimVector, blockDimVector>>>(Pb, P, N);
  cudaThreadSynchronize();

  /*
  double *temp;
  cudaMalloc((void**)&temp, N*sizeof(double));

  thrust::device_ptr<double> dev_ptr(temp);

  for (int i=0; i<N;i++){
    kernelForwardMap<<<gridDimVector, blockDimVector>>>(temp, L, y, i,N);
    cudaThreadSynchronize();

    double result = thrust::reduce(dev_ptr, dev_ptr + i);

    kernelForwardSolve<<<1,1>>>(Pb,y,result,i,N);
    cudaThreadSynchronize();
  }
  cudaFree(temp);*/

  kernelSequentialHelpAfter<<<1,1>>>(A,L,U,P,Pb,b,x,y,N);
  cudaThreadSynchronize();
}


std::vector<std::tuple<int, int>> matching(double* graph, int n, double err){

  int matrix_size = n;
  std::vector<std::tuple<int, int>> M = std::vector<std::tuple<int, int>>();

  vector<int> rc_map = vector<int>();
  for (int i=0; i<n; i++){
    rc_map.push_back(i);
  }

  double *host_x;
  cudaMallocHost((void **) &host_x, n*sizeof(double));

  //double* host_x = (double *) calloc(n, sizeof(double));
  while (M.size() < n/2){
      //kernelSequentialSolve<<<1,1>>>(A,L,U,P,Pb,b,x,y,matrix_size);
      //cudaThreadSynchronize();
      solve(matrix_size);

      int row_j = -1;
      cudaMemcpy(host_x, x, matrix_size*sizeof(double), cudaMemcpyDeviceToHost);

      for (int row=0; row< matrix_size; row++){
        int true_first_col = rc_map[0];
        int true_row = rc_map[row];
        if (abs(host_x[row]) > err && graph[true_row*n + true_first_col]!=0){
          row_j = row;
          break;
        }
      }

      if (row_j == -1){
        printf("row_j is -1\n");
        printf("matching size = %d\n", M.size());
        return std::vector<std::tuple<int, int>>();
      }
      /*
      printf("\nnew edge: (%d, %d)\n", 0, row_j);
      printf("host_x[row]: %.3e\n", host_x[row_j]);
      printf("host abs condition: %s\n", fabs(host_x[row_j]) > err ? "true" : "false");
      printf("graph edge: %.3e\n", graph[rc_map[row_j]][rc_map[0]]);
  */
      M.push_back(std::make_tuple(rc_map[0], rc_map[row_j]));
      int row_counter = 0;

      for (int i=0; i<matrix_size; i++){
        if (i != 0 && i != row_j){

          dim3 blockDim(1024,1,1);
          dim3 gridDim((matrix_size + blockDim.x - 1)/blockDim.x,1,1);

          kernelresizeA<<<gridDim, blockDim>>>(A, row_counter, i, row_j, matrix_size);
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

std::vector<std::tuple<int, int>> setup(double *cudaGraph, vector<vector<double>> host_graph, int N, double err){
    
    printf("error bound is: %.3e\n\n", err);
    dim3 blockDim(1024,1,1);
    dim3 gridDim((N + blockDim.x - 1)/blockDim.x,1,1);

    cudaMalloc((void**)&A, N*N*sizeof(double));
    cudaMalloc((void**)&P, sizeof(double)*N*N);
    cudaMalloc((void**)&L, sizeof(double)*N*N);
    cudaMalloc((void**)&U, sizeof(double)*N*N);
    cudaMalloc((void**)&Pb, sizeof(double)*N);
    cudaMalloc((void**)&b, N*sizeof(double));
    cudaMalloc((void**)&x, N*sizeof(double));
    cudaMalloc((void**)&y, N*sizeof(double));

    kernelSetb<<<gridDim, blockDim>>>(b, N);
    cudaThreadSynchronize();

    cudaMemcpy(A, cudaGraph, N*N*sizeof(double), cudaMemcpyHostToDevice);
    return matching(cudaGraph, N, err);  
}

/*

__global__ void  kernelSequentialSolve(double *A, double *L, double *U, double *P,
  double *Pb, double *b, double *x, double *y, int N){


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
    if (P[i*N] == 1.0){
      Pb[i] = 1.0;
    } else {
      Pb[i] = 0.0;
    }
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


*/



