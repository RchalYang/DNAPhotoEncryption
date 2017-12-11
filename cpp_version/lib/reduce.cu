#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdint>
#include <iostream>
#include "reduce.h"
#include "dna.h"

__global__ void reduceKernel(int *d_in, int* d_out, int totalsize) {

	int myID = threadIdx.x + blockIdx.x * blockDim.x;
	int tid = threadIdx.x;

	extern __shared__ float sdata[];

	if (myID < totalsize) {
        sdata[tid] = d_in[myID];
	}
	else {
		sdata[tid] = 0;
	}

	__syncthreads();

	if (myID >= totalsize) {
		if (tid == 0) {
			d_out[blockIdx.x] = 0;
		}
		return;
	}
    
	//divide threads into two parts according to threadID, and add the right part to the left one, lead to reducing half elements, called an iteration; iterate until left only one element
	for (unsigned int s = blockDim.x / 2; s>0; s >>= 1) {
		if (tid < s) {
			sdata[tid] = sdata[tid + s] + sdata[tid];
		}
		__syncthreads(); //ensure all adds at one iteration are done
	}

	if (tid == 0) {
		d_out[blockIdx.x] = sdata[0];
	}
}

__global__ void initialKernel(uint8_t* d_input, int* d_output, int totalsize, int toReduce){

	int myID = threadIdx.x + blockIdx.x * blockDim.x;

	if (myID < totalsize) {
        if (d_input[myID] == toReduce){
		    d_output[myID] = 1;
        }  
        else{
            d_output[myID] = 0;
        }
	}

};

void reduction(const uint8_t* const h_in, int* result, int input_size, int toReduce) {
	int blocksize = 64;

	uint8_t* d_temp;
	cudaMalloc(&d_temp, input_size * sizeof(uint8_t));
	cudaMemcpy(d_temp, h_in, input_size * sizeof(uint8_t), cudaMemcpyHostToDevice);

    int* d_raw;
    cudaMalloc(&d_raw, input_size * sizeof(int));
	
    initialKernel<<<roundup(input_size, blocksize),blocksize>>>(d_temp, d_raw, input_size, toReduce);

    cudaFree(d_temp);

	const int shared_mem_size = sizeof(int)*blocksize;
	int currentSize = input_size;
	int tsize = input_size;
	while (1) {
		int* d_current_out;
		tsize = currentSize;
		currentSize = roundup(tsize, blocksize);
		// std::cout << currentSize << std::endl;
		cudaMalloc(&d_current_out, currentSize * sizeof(int));

		reduceKernel <<< currentSize, blocksize, shared_mem_size >>> (d_raw, d_current_out, tsize);
		cudaFree(d_raw);

		d_raw = d_current_out;

		// int* temp_out = (int*)malloc(sizeof(int)*currentSize);
		// cudaMemcpy(temp_out, d_raw, currentSize * sizeof(int), cudaMemcpyDeviceToHost);
		// for (int p = 0; p < currentSize; p++) {
		// 	std::cout << temp_out[p]<<" ";
		// }
		// free(temp_out);
		// std::cout << std::endl;

		if (currentSize == 1)
		{
			break;
		}
	}

	cudaMemcpy(result, d_raw, sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(d_raw);
}
