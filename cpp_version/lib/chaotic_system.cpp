#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp> 
#include <iostream>
#include <vector>

float logistic_func(float u, float x){
    return u*x*(1-x);
}

void logistic_map(float* map, float u, float x,int size){
    map[0] = x;
    for (int i=1; i<size; i++){
        x = logistic_func(u,x);
        map[i] = x;
    }
}


float spatiotemporal_func(int idx, float* map_previous, double u, float e){
    int ip1 = (idx+1)%3;
    int im1 = (idx+2)%3;
    return (1-e)*logistic_func(u, map_previous[idx]) + e / 2 * (logistic_func(u, map_previous[ip1])+logistic_func(u, map_previous[im1])) ;
}


void spatiotemporal_system(float* output_map, float* map_zero, float e, double u,int size){
    output_map[0] = spatiotemporal_func(0, map_zero, u, e);
    output_map[1] = spatiotemporal_func(1, map_zero, u, e);
    output_map[2] = spatiotemporal_func(2, map_zero, u, e);

    for (int i=1; i < size; i++){
        float* temp_map_ptr = output_map + (i-1)*3;
        output_map[i*3] = spatiotemporal_func(0, temp_map_ptr, u, e);
        output_map[i*3+1] = spatiotemporal_func(1, temp_map_ptr, u, e);
        output_map[i*3+2] = spatiotemporal_func(2, temp_map_ptr, u, e);
    }
    
}
