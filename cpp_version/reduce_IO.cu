#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp> 
#include <iostream>
#include <chrono>
#include "./lib/chaotic_system.h"
#include "./lib/config.h"
#include "./lib/transition.h"
#include "./lib_io_reduce/dna_io_reduce.h"
#include "./lib_io_reduce/reduce_io_reduce.h"
#include "./lib_io_reduce/rotate_io_reduce.h"
#include <cmath>
#include <string>
#include <algorithm>

#define EPS 0.00000000000000001

bool isZero(double a){
    return ( fabs(a) < 1e-18);
}


int compare(const void* a, const void* b){
    return ((RC*)a)->value > ((RC*)b)->value;
}


int test( cv::Mat image ) {

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    
    int cols = image.cols;
    int rows = image.rows;
    int LargeSize = 4*cols*rows;
    int SmallSize = cols*rows;

    std::cout<<"Resolution:"<<SmallSize<<std::endl;

    uint8_t* h_image = (uint8_t*) malloc(sizeof(uint8_t)*SmallSize);

    uint8_t* image_data;
    for (int row = 0; row<rows; row++){
        image_data = image.ptr<uint8_t>(row);
        for (int col = 0; col < cols; col++){
            h_image[row*cols+col] = image_data[col];
        }
    }

    float* logi_map = (float*) malloc(sizeof(float)*LargeSize);
    logistic_map(logi_map,u0_l,x0_l,LargeSize);

    float* d_logi_map;
    cudaMalloc(&d_logi_map,sizeof(float)*LargeSize);
    cudaMemcpy(d_logi_map, logi_map, LargeSize * sizeof(float), cudaMemcpyHostToDevice);

    uint8_t* d_image;
    cudaMalloc(&d_image,sizeof(uint8_t)*SmallSize);
    cudaMemcpy(d_image, h_image, SmallSize * sizeof(uint8_t), cudaMemcpyHostToDevice);

    uint8_t* logi_map_dna = image_to_dna_matrix_float(d_logi_map, LargeSize);

    uint8_t DNArules[] = {  1,1,0,0,3,3,2,2,
                            0,3,1,2,1,2,0,3,
                            3,0,2,1,2,1,3,0,
                            2,2,3,3,0,0,1,1 };


    
    d_image = image_to_dna_matrix_uint8(d_image, d_logi_map, DNArules, SmallSize, 4*cols);

    free(logi_map);
    cudaFree(d_logi_map);

    uint8_t* d_dna_map = dna_add(logi_map_dna,d_image,4*rows*cols);

    int sum[4];
    for (int i=0;i<4;i++){
        reduction(d_dna_map,&sum[i],4*rows*cols,i);
        // reduction()
        std::cout<<sum[i]<<std::endl;
    }

    double H=0;
    for (int i=0;i<4;i++){
        double p = sum[i]/(1.0*rows*cols*4);
        H += log2(1.0/p)*p;
    }

    double x0_h = H - floorf(H);
    // std::cout<<x0_h<<std::endl;
    
    if (isZero(x0_h)){
        x0_h = x0_l;
    }

    double x20_h = x0_h;
    double x10_l = x0_l;
    for (int i=0;i<10;i++){
        x20_h = logistic_func(u0_l,x20_h);        
        x20_h = logistic_func(u0_l,x20_h);

        x10_l = logistic_func(u0_l,x10_l);
    }

    // std::cout<<x20_h<<std::endl;

    std::string h_bin = dec2bin(x20_h);
    double h_dec = bin2dec(h_bin);

    // std::cout<<h_dec<<std::endl;

    // std::cout<<x10_l<<std::endl;

    std::string tempstr = dec2bin(x10_l);

    // std::cout<<tempstr<<std::endl;

    std::string h_cipher = xorstr(h_bin, tempstr);
    
    std::cout<<h_cipher<<std::endl;

    double u_s = 3.75 + 0.25 * h_dec;

    int max_m_4n = std::max(rows,4*cols);
    float* sc_map = (float*)malloc(sizeof(float)*max_m_4n*3);

    float x0_s[3] = {0.364,0.785,0.293};
    spatiotemporal_system(sc_map, x0_s, epsilon, u_s, max_m_4n);
    // spatiotemporal_system

    RC* R = (RC*)malloc(sizeof(RC)*rows);
    RC* C = (RC*)malloc(sizeof(RC)*4*cols);
    
    for (int i=0; i<rows; i++){
        R[i].value=sc_map[3*i];
        R[i].idx = i;
    }

    for (int i=0; i<4*cols; i++){
        C[i].value=sc_map[3*i+2];
        C[i].idx = i;
    }

    qsort(R,rows,sizeof(RC),compare);
    qsort(C,4*cols,sizeof(RC),compare);

    rotate(d_dna_map, rows, cols, R, C);

    uint8_t* dna_map = dna_matrix_to_image(d_dna_map, cols*rows);

    
    for (int row = 0; row<rows; row++){
        image_data = image.ptr<uint8_t>(row);
        for (int col = 0; col < cols; col++){
            image_data[col] = dna_map[row*cols+col];
        }
    }


    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double ttrack1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2-t1).count();

    // cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
    // cv::imshow( "Display window", image );
    cv::imwrite("result_io_reduce.jpg", image);
    
    std::cout<<"Total Time Consumption:"<<ttrack1<<std::endl;
    free(dna_map);

    cv::waitKey(0);
    return 0;
}

int main(int argc, char** argv){
    cv::Mat image;
    image = cv::imread("lena30.jpg", CV_LOAD_IMAGE_GRAYSCALE);
    
    if(! image.data ) {
        std::cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }

    test()
    

}