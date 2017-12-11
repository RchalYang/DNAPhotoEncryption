#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp> 
#include <iostream>
#include <chrono>
#include "./lib/chaotic_system.h"
#include "./lib/config.h"
#include "./lib/transition.h"
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

int roundup(int a, int b){
    return ceil((1.0*a)/ (1.0* b));
}

int main( int argc, char** argv ) {
  
    double totaltime1 = 0;

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    
    cv::Mat image;
    image = cv::imread("lena30.jpg", CV_LOAD_IMAGE_GRAYSCALE);
    
    if(! image.data ) {
        std::cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }

    int cols = image.cols;
    int rows = image.rows;

    int LargeSize = rows*cols*4;
    int SmallSize = rows*cols;

    uint8_t* h_image = (uint8_t*) malloc(sizeof(uint8_t)*rows*cols);

    uint8_t* image_data;
    for (int row = 0; row<rows; row++){
        image_data = image.ptr<uint8_t>(row);
        for (int col = 0; col < cols; col++){
            h_image[row*cols+col] = image_data[col];
        }
    }

    float* logi_map = (float*) malloc(sizeof(float)*LargeSize);

    logistic_map(logi_map,u0_l,x0_l,LargeSize);

    uint8_t* logi_map_dna = (uint8_t*) malloc(sizeof(uint8_t)*LargeSize);

    for (int myIdx=0;myIdx<LargeSize;myIdx++){
        int temp = (int(logi_map[myIdx]*256.0))%4;
        switch (temp){
            case 0:
                logi_map_dna[myIdx] = 1;
                break;
            case 1:
                logi_map_dna[myIdx] = 0;
                break;
            case 2:
                logi_map_dna[myIdx] = 3;
                break;
            case 3:
                logi_map_dna[myIdx] = 2;
        }
    }

    // std::cout<<"axiba"<<std::endl;

    uint8_t DNArules[] = { 1,1,0,0,3,3,2,2,
                            0,3,1,2,1,2,0,3,
                            3,0,2,1,2,1,3,0,
                            2,2,3,3,0,0,1,1 };

    // h_image = image_to_dna_matrix_uint8(h_image, logi_map, DNArules, rows*cols, 4*cols);

    uint8_t* tOutput = (uint8_t*) malloc(sizeof(uint8_t)*LargeSize);

    for (int myIdx=0;myIdx<SmallSize;myIdx++){

        int inputIdx = myIdx/4;
        uint8_t temp_input = h_image[myIdx];
        uint8_t temp_flat = (int)(logi_map[roundup(myIdx,4*cols)*4*cols]*256)%8;
        
        uint8_t tempIdx = myIdx - 4 * inputIdx + 1;
        temp_input = temp_input / (1<<(2*(4-tempIdx))) % 4;

        tOutput[myIdx] = DNArules[temp_input*8+temp_flat];
    }

    std::cout<<"axiba"<<std::endl;

    free(logi_map);
    free(h_image);
    h_image = tOutput;

    
    for (int myIdx=0;myIdx<LargeSize;myIdx++){
        h_image[myIdx] = (h_image[myIdx]+logi_map_dna[myIdx])%4;
    }

    std::cout<<"axiba"<<std::endl;
    
    uint8_t* dna_map = h_image;

    free(logi_map_dna);


    int sum[4];
    sum[0]=0;
    sum[1]=0;
    sum[2]=0;
    sum[3]=0;

    for (int myIdx=0;myIdx<LargeSize;myIdx++){
        sum[dna_map[myIdx]]++;
    }

    for (int i=0;i<4;i++){
        std::cout<<sum[i]<<std::endl;
    }

    double H=0;
    for (int i=0;i<4;i++){
        double p = sum[i]/(1.0*rows*cols*4);
        H += log2(1.0/p)*p;
    }

    double x0_h = H - floorf(H);
    std::cout<<x0_h<<std::endl;
    
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

    std::string h_bin = dec2bin(x20_h);
    double h_dec = bin2dec(h_bin);

    std::string tempstr = dec2bin(x10_l);

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

    // rotate(dna_map, rows, cols, R, C);
    tOutput = (uint8_t*) malloc(sizeof(uint8_t)*LargeSize);

    for (int myIdx=0;myIdx<LargeSize;myIdx++){
        int tRow = myIdx/(4*cols);
        int tCol = myIdx%(4*cols);
        tOutput[myIdx] = dna_map[(R[tRow].idx)*cols+C[tCol].idx];
        // h_image[myIdx] = (h_image[myIdx]+logi_map_dna[myIdx])%4;
    }    

    free(dna_map);
    dna_map = tOutput;

    tOutput = (uint8_t*) malloc(sizeof(uint8_t)*SmallSize);

    // dna_map = dna_matrix_to_image(dna_map, cols*rows);
    for (int myIdx=0;myIdx<SmallSize;myIdx++){
        uint8_t temp_base = 64;
        uint8_t temp_data = 0;
        for (int i=0;i<4;i++,temp_base/=4){
            temp_data += ((5-dna_map[4*myIdx + i])%4)*temp_base;
        }
        tOutput[myIdx] = temp_data;
    }

    free(dna_map);
    dna_map = tOutput;

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
    cv::imwrite("result_Sequencial.jpg", image);
    
    std::cout<<"Total Time Consumption:"<<ttrack1<<std::endl;
    free(dna_map);

    cv::waitKey(0);
    return 0;
}