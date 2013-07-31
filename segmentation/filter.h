#include "gammaTone.h"

class filter
{
	int fs;
	int fLength, nOrder, nFFT;

	float beta;

	float *rValue;
	float *iValue;

	float *inputR;
	float *inputI;

	void hammingPara(float delta, float transBw, int &fLength);

	void hammingBandPass(float wc, float wn);
	void hammingDiffOneBandPass(float wc, float wn);
	void hammingDiffTwoBandPass(float wc, float wn);

	void kaiserPara(float delta, float transBw);

	void kaiserLowPass(float wn);
	void diffKaiserLowPass(float wn);

public:
	filter(float ws, float wn, float ts, int fs, int type);

	~filter();

	void filtering(float *input, float *output, int sigLength);
};