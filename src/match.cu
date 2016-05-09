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

#define THREADS_PB 8

__device__ double *A;
__device__ double *L;
__device__ double *U;
__device__ double *P;
__device__ double *Pb;
__device__ double *b;
__device__ double *x;
__device__ double *y;
__device__ int *i;
__device__ double *temp;


__global__ void kernelRowSwap(double *M, int start, int end, int k, int i, int N){
    
   
    int col = blockDim.x*blockIdx.x + threadIdx.x;

    if (col < start || col >= end)
        return;

    double temp = M[k*N + col];
    M[k*N + col] = M[i*N + col];
    M[i*N + col] = temp;
}

__global__ void kernelSetL(double *L, double *U, int k, int N){
  int j = blockDim.x*blockIdx.x + threadIdx.x;
  if (j >= N || j<k+1)
        return;

  L[j*N+k] = U[j*N + k]/U[k*N+k];
}

__global__ void kernelSetLNew(double *U, int k, int N){
  int j = blockDim.x*blockIdx.x + threadIdx.x;
  if (j >= N || j<k+1)
        return;

  U[j*N+k] = U[j*N + k]/U[k*N+k];
}

__global__ void kernelSetU(double *L, double *U, int k, int N){
    int j = blockDim.x*blockIdx.x + threadIdx.x;
    int col = blockIdx.y * blockDim.y + threadIdx.y;

    if (j >= N || j<k+1 || col >= N || col <k)
        return;

    U[j*N + col] -= L[j*N+k]*U[k*N+col];
}

__global__ void kernelSetUNew(double *U, int k, int N){
    int j = blockDim.x*blockIdx.x + threadIdx.x;
    int col = blockIdx.y * blockDim.y + threadIdx.y;

    if (j >= N || j<k+1 || col >= N || col <k+1)
        return;

    
    double temp = U[j*N+k]*U[k*N+col];
    U[j*N + col] -= temp;
}


__global__ void kernelInitL(double *L, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    L[i*N+j] = (i==j) ? 1.0 : 0.0;
}

__global__ void kernelInitP(double *P, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N)
        return;

    P[i] = (i==0) ? 1.0 : 0.0;
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

__global__ void kernelFillTemp(double *temp, double *U, int k, int N){

  int row = blockDim.x*blockIdx.x + threadIdx.x;
  if (row >= N)
    return;

  if (row >= k){
    temp[row] = abs(U[row*N+k]);
  } else {
    temp[row] = 0.0;
  }
}

__global__ void kernelSwapP(double *P, int k, int i, int N){

  double temp = P[k];
  P[k] = P[i];
  P[i] = temp;
}

__global__ void kernelGetUfromA(double *U, double*A, int N){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= N || j >= N || i > j)
    return;

  U[i*N+j] = A[i*N+j];
}

__global__ void kernelGetLfromU(double *L, double*U, int N){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= N || j >= N || i < j)
    return;

  if (i==j){
    L[i*N+j] = 1.0;
  } else {
    L[i*N+j] = U[i*N+j];
  }
}

void forward_solve(int N){

  dim3 blockDimVector(1024,1,1);
  dim3 gridDimVector((N + blockDimVector.x - 1)/blockDimVector.x,1,1);

  thrust::device_vector<double> PB_vec(Pb, Pb+N);
  thrust::device_vector<double> y_vec(y, y+N);
  thrust::device_vector<double> U_vec(U, U+(N*N));
  thrust::device_vector<double> d_vec(temp, temp+N);

  for (int i = 0; i<N; i++){
    for (int j=0; j<i; j++){
        d_vec[j] = U_vec[i*N+j]*y_vec[j];
    }
    double sum = thrust::reduce(d_vec.begin(), d_vec.begin()+i);

    y_vec[i] = PB_vec[i] - sum;
    printf("fsolve val: %.3e\n", y_vec[i]);
  }
  printf("\n");
}

__global__ void kernelSequentialHelpAfter(double *A, double *U, 
  double *Pb, double *x, double *y, int N){

  for (int i = 0; i < N; i++){
    double rhs = Pb[i];
    for (int j=0; j<i; j++){
      rhs -= U[i*N+j]*y[j];
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
  kernelInitP<<<gridDimVector, blockDimVector>>>(Pb,N);
  //cudaDeviceSynchronize();

  thrust::device_vector<double> d_vec(N);
  double* pd_vec = thrust::raw_pointer_cast(d_vec.data());
  thrust::device_vector<double>::iterator iter;

  for (int k=0; k<N-1; k++){

    kernelFillTemp<<<gridDimVector, blockDimVector>>>(pd_vec, U, k, N);
    //cudaDeviceSynchronize();
    iter = thrust::max_element(d_vec.begin(), d_vec.end());
    unsigned int position = iter - d_vec.begin();
    kernelRowSwap<<<gridDimVector, blockDimVector>>>(U, 0, N, k,position, N);
    kernelSwapP<<<1,1>>>(Pb,k,position,N);
    //cudaDeviceSynchronize();

    kernelSetLNew<<<gridDimVector, blockDimVector>>>(U,k,N);
    //cudaDeviceSynchronize();
    kernelSetUNew<<<gridDimArray, blockDimArray>>>(U,k,N);
    //cudaDeviceSynchronize();
  }

  kernelSequentialHelpAfter<<<1,1>>>(A,U,Pb,x,y,N);
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
      //cudaDeviceSynchronize();
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
          //cudaDeviceSynchronize();
          row_counter++;
        }
      }
      matrix_size -= 2;
      rc_map.erase(rc_map.begin() + row_j);
      rc_map.erase(rc_map.begin());
  }
  return M;
}

double* init_cudaGraph(int N){
  double *cudaGraph;
  cudaMallocHost((void**)&cudaGraph, N*N*sizeof(double));
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      cudaGraph[i*N+j] = 0.0;
    }
  }
  return cudaGraph;
}

std::vector<std::tuple<int, int>> setup(double *cudaGraph, vector<vector<double>> host_graph, int N, double err){
    
    printf("error bound is: %.3e\n\n", err);
    dim3 blockDim(1024,1,1);
    dim3 gridDim((N + blockDim.x - 1)/blockDim.x,1,1);

    cudaMalloc((void**)&A, N*N*sizeof(double));
    cudaMalloc((void**)&U, sizeof(double)*N*N);
    cudaMalloc((void**)&Pb, sizeof(double)*N);
    cudaMalloc((void**)&x, N*sizeof(double));
    cudaMalloc((void**)&y, N*sizeof(double));
    cudaMalloc((void**)&i, sizeof(int));
    cudaMemcpy(A, cudaGraph, N*N*sizeof(double), cudaMemcpyHostToDevice);

    kernelInitP<<<gridDim, blockDim>>>(Pb, N);
    //cudaDeviceSynchronize();

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

/*
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
*/


