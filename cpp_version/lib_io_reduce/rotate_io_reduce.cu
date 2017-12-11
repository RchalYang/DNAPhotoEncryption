#include "rotate_io_reduce.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "../lib/config.h"
#include <cstdint>
#include <iostream>

__global__ void rotate_kernel(uint8_t* d_input, uint8_t* d_output, int rows, int cols, RC* d_R, RC* d_C){

    int myID = blockIdx.x*blockDim.x + threadIdx.x;
    int rowID = myID/(4*cols);
    int colID = myID%(4*cols);

    if (myID >= rows*cols*4 ){
        return;
    }
    
    int R_target = d_R[rowID].idx;
    int C_target = d_C[colID].idx;
    d_output[myID] = d_input[R_target*(4*cols)+C_target];
}

void rotate(uint8_t* d_dna_map, int rows, int cols, RC* R, RC* C){
    int blockSize = 64;
    
    uint8_t* d_out;
    cudaMalloc(&d_out,sizeof(uint8_t)*4*cols*rows);
    
    RC* d_R;
    RC* d_C;
    cudaMalloc(&d_R,sizeof(RC)*rows);
    cudaMalloc(&d_C,sizeof(RC)*4*cols);
    cudaMemcpy(d_R, R, sizeof(RC)*rows, cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, sizeof(RC)*4*cols, cudaMemcpyHostToDevice);

    rotate_kernel<<<ceil((4.0*cols*rows)/ (1.0* blockSize)),blockSize >>> (d_dna_map,d_out,rows,cols,d_R,d_C);

    cudaFree(d_R);
    cudaFree(d_C);
    cudaFree(d_dna_map);
    d_dna_map = d_out;
}