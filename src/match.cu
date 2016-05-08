#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

#define THREADS_PB 256


__global__ void kernelRowSwap(double *M, int start, int end, int k, int i, int N){
    int threadIndex = blockDim.x*blockIdx.x + threadIdx.x;

    if (threadIndex >= (start - end))
        return;

    double temp = M[k*N + start + threadIndex];
    M[k*N + start + threadIndex] = M[i*N + start + threadIndex];
    M[i*N + start + threadIndex] = temp;
}



__global__ void kernelSetLURow(double *L, double *U, int k, int N){
    int j = blockDim.x*blockIdx.x + threadIdx.x;

    //ERROR !!!! check this
    int m = 0;

    if (j >= N || j<k+1)
        return;

    L[j*N + k] = U [j*N + k]/U[k*N+k];
    for (int col=k; col<m; col++){
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

__global__ void copyAtoU(double *A, double *U, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N || j >= N)
        return;

    U[i*N+j] = .1*A[i*N+j];
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


void setup(){
    int N = 100;
    dim3 blockDim(1024,1,1);
    dim3 gridDim((N + blockDim.x - 1)/blockDim.x,1,1);


    int start = 1;
    int end = 1;
    int k = 1;
    int i = 1;
    double *M;
    cudaMalloc((void**)&M, sizeof(double)*N);


    kernelRowSwap<<<gridDim, blockDim>>>(M,start,end,k,i,N);
    cudaThreadSynchronize();

    double *L;
    double *U;
    cudaMalloc((void**)&L, sizeof(double)*N);
    cudaMalloc((void**)&U, sizeof(double)*N);


    kernelSetLURow<<<gridDim, blockDim>>>(L, U, k,  N);

    double *I;
    cudaMalloc((void**)&I, sizeof(double)*N);

    kernelMakeIdentity<<<gridDim, blockDim>>>(I, N);


    double *A;
    cudaMalloc((void**)&A, sizeof(double)*N);

    copyAtoU<<<gridDim,blockDim>>>(A, U,  N);

    int row_j = 0;
    int row_counter = 0;
    kernelresizeA<<<gridDim, blockDim>>>(A, row_counter, i, row_j,  N);

    double *Pb;
    cudaMalloc((void**)&Pb, sizeof(double)*N);

    double *P;
    cudaMalloc((void**)&P, sizeof(double)*N);


    kernelSetPb<<<gridDim, blockDim>>>(Pb,P, N);
}


