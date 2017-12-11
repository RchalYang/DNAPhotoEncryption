#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdint>

uint8_t* dna_add(uint8_t* a, uint8_t* b, int size);

uint8_t* image_to_dna_matrix_uint8(uint8_t*,float*,uint8_t*, int, int);

uint8_t* image_to_dna_matrix_float(float*, int);

uint8_t* dna_matrix_to_image(uint8_t*, int);

int roundup(int, int);

//A=1  C=0 G=3 T=2