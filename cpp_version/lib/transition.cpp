#include "transition.h"
#include <string>
#include <iostream>

double bin2dec(std::string source){

    double result = 0;
    double temp = 1.0/2;
    for (int i = 0; i<64;i++){
        result = result + temp*(source[i]-48);
        temp /= 2;
    }
    return result;
}

std::string dec2bin(double source){
    std::string result("");

    for (int i=0;i<64;i++){
        source *= 2;
        // result <<= 1;
        if ((source-1) > 1e-17 ){
            source -= 1;
            result.append("1"); 
        }else{
            result.append("0");
        }
    }

    return result;

}


std::string xorstr(std::string left, std::string right){
    std::string result("");
    for (int i=0;i<64;i++){
        char temp = (left[i]-48) ^ (right[i]-48);
        // std::string temp_S(temp);
        result.push_back(temp+48);
        // std::cout<<i<<std::endl;
    }
    return result;
}

// function result=my_dec2bin(dec,pre)
// result = '0.';
// for i = 1:pre
//     dec = dec*2;
//     if dec>1
//         result = strcat(result,'1');
//         dec = dec-1;
//     else
//         result = strcat(result,'0');
//     end
// end

// function result=my_bin2dec(bin)
// result = 0;
// for i = 1:64
//     result = result + str2num(bin(i+2))*2^(-i);
// end