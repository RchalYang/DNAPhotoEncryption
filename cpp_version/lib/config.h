#ifndef CONFIG
#define CONFIG
const float u0_l = 3.785;
const float x0_l = 0.437;

// float x0_s[3] = {0.364,0.785,0.293};

const float epsilon = 0.2582;

struct _RC{
    float value;
    int idx;
};

typedef struct _RC RC;

#endif