#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI (3.1415926535897932384626433832795)

inline float sigmoid(float x){ return( 2/(1+exp(-2*x)) -1 ); }

inline float hamming(float x){ return( 0.54 + 0.46*cos(x*PI) ); }

inline float hammingDiffOne(float x){ return( -0.46*PI*sin(x*PI) ); }

inline float hammingDiffTwo(float x){ return( -0.46*PI*PI*cos(x*PI) ); }

float sinc(float tim, float wn);
float sincDiffOne(float tim, float wn);
float sincDiffTwo(float tim, float wn);

float bessi0(float x);

float bessi1(float x);

int maxPos(float *input, int d1, int d2);

int minPos(float *input, int d1, int d2);

void fft(float *inputR, float *inputI, int N, float direct);

float zeroCross(float *input, int sLen);

int findLocalPeak(float *input, int p1, int p2, float theta);

int findLocalValley(float *input, int p1, int p2, float theta);