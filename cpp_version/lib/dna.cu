#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdint>
#include <iostream>
#include "dna.h"
using namespace std;

__device__ int roundup_d(int a, int b){
    return ceil((1.0*a)/ (1.0* b));
}

int roundup(int a, int b){
    return ceil((1.0*a)/ (1.0* b));
}


__global__ void dna_add_kernel(uint8_t* d_in_left, uint8_t* d_in_right, int size){
    int myIdx = blockIdx.x*blockDim.x + threadIdx.x;
    if (myIdx < size){
        d_in_left[myIdx] = (d_in_right[myIdx] + d_in_left[myIdx])%4;
    }
}


uint8_t* dna_add(uint8_t* a, uint8_t* b, int size){
    uint8_t* d_left;
    uint8_t* d_right;
    cudaMalloc(&d_left,sizeof(uint8_t)*size);
    cudaMalloc(&d_right,sizeof(uint8_t)*size);
    cudaMemcpy(d_left, a, size*sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_right, b, size*sizeof(uint8_t), cudaMemcpyHostToDevice);

    int blockSize = 64;

    dna_add_kernel<<<roundup(size,blockSize),blockSize>>>(d_left,d_right,size);

    // uint8_t* output = (uint8_t*) malloc(sizeof(uint8_t)*size);
    cudaMemcpy(a, d_left, size*sizeof(uint8_t), cudaMemcpyDeviceToHost);

    cudaFree(d_left);
    cudaFree(d_right);
    free(b);
    return a;
}


__global__ void dna_encode_kernel_lmatrix(float* d_input, uint8_t* d_output, int size){
    int myIdx = blockIdx.x*blockDim.x + threadIdx.x;

    if (myIdx >= size)
        return;

    float temp_data = d_input[myIdx];
    int temp = (int(temp_data*256.0))%4;
    switch (temp){
        case 0:
            d_output[myIdx] = 1;
            break;
        case 1:
            d_output[myIdx] = 0;
            break;
        case 2:
            d_output[myIdx] = 3;
            break;
        case 3:
            d_output[myIdx] = 2;
    }

}


__global__ void dna_encode_kernel_image(uint8_t* d_input, float* reference, uint8_t* dna_rules, uint8_t* d_output, int size, int cols){
    int myIdx = blockIdx.x*blockDim.x + threadIdx.x;
	int inputIdx = myIdx/4;

    if (inputIdx >= size)
        return;

    uint8_t temp_input = d_input[inputIdx];
    uint8_t temp_flat = (int)(reference[roundup_d(myIdx,cols)*cols]*256)%8;
    
    // __synthreads();
	uint8_t tempIdx = myIdx - 4 * inputIdx + 1;
	temp_input = temp_input / (1<<(2*(4-tempIdx))) % 4;

    d_output[myIdx] = dna_rules[temp_input*8+temp_flat];
}


uint8_t* image_to_dna_matrix_uint8(uint8_t* input, float* reference, uint8_t* h_rules, int input_size, int cols){

	int blockSize = 64;
    int ruleSize = 32;

    uint8_t* d_image;
    uint8_t* d_matrix;
	float* d_reference;
    uint8_t* d_rules;
    cudaMalloc(&d_image,sizeof(uint8_t)*input_size);
    cudaMalloc(&d_matrix,sizeof(uint8_t)*4*input_size);
    cudaMalloc(&d_reference,sizeof(float)*4*input_size);
    cudaMalloc(&d_rules,sizeof(uint8_t)*ruleSize);
    cudaMemcpy(d_image, input, input_size * sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_reference, reference, 4 * input_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rules, h_rules, ruleSize*sizeof(uint8_t), cudaMemcpyHostToDevice);
    
	dna_encode_kernel_image <<< roundup(4*input_size,blockSize), blockSize >>> (d_image, d_reference, d_rules, d_matrix,input_size, cols);

    uint8_t* output = (uint8_t*) malloc( sizeof(uint8_t)*4*input_size);
    cudaMemcpy(output, d_matrix, sizeof(uint8_t)*4*input_size, cudaMemcpyDeviceToHost);
    cudaFree(d_image);
    cudaFree(d_matrix);
    free(input);

    return output;
};

uint8_t* image_to_dna_matrix_float(float* input, int input_size){
	
    int blockSize = 64;

    float* d_image;
    uint8_t* d_matrix;
    cudaMalloc(&d_image,sizeof(float)*input_size);
    cudaMalloc(&d_matrix,sizeof(uint8_t)*input_size);
    cudaMemcpy(d_image, input, input_size * sizeof(float), cudaMemcpyHostToDevice);

    cout<<roundup(input_size,blockSize)<<" "<<blockSize<<endl;

    dna_encode_kernel_lmatrix <<< roundup(input_size,blockSize), blockSize >>> (d_image,d_matrix,input_size);

    uint8_t* output = (uint8_t*) malloc( sizeof(uint8_t)*input_size);
    cudaMemcpy(output, d_matrix, sizeof(uint8_t)*input_size, cudaMemcpyDeviceToHost);
    cudaFree(d_image);
    cudaFree(d_matrix);

    return output;
};

__global__ void dna_matrix_decode_kernel(uint8_t* d_input, uint8_t* d_output, int size){
    int myIdx = blockIdx.x*blockDim.x + threadIdx.x;

    if (myIdx >= size)
        return;

    uint8_t temp_base = 64;
    uint8_t temp_data = 0;
    for (int i=0;i<4;i++,temp_base/=4){

        temp_data += ((5-d_input[4*myIdx + i])%4)*temp_base;
    }

    d_output[myIdx] = temp_data;
    
}

uint8_t* dna_matrix_to_image(uint8_t* input_matrix, int size){
    int blockSize = 64;

    uint8_t* d_image;
    uint8_t* d_matrix;
    cudaMalloc(&d_image,sizeof(uint8_t)*size);
    cudaMalloc(&d_matrix,sizeof(uint8_t)*4*size);
    cudaMemcpy(d_matrix, input_matrix, 4 * size * sizeof(uint8_t), cudaMemcpyHostToDevice);

    dna_matrix_decode_kernel <<< roundup(size,blockSize), blockSize >>> (d_matrix,d_image,size);

    uint8_t* output = (uint8_t*) malloc( sizeof(uint8_t)*size);
    cudaMemcpy(output, d_image, sizeof(uint8_t)*size, cudaMemcpyDeviceToHost);
    cudaFree(d_image);
    cudaFree(d_matrix);

    free(input_matrix);

    return output;
}