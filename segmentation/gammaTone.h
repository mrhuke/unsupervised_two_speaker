#include "tool.h"
#include "common.h"

#define BW_CORRECTION       1.019      			/* ERB bandwidth correction 4th order */

struct gammaTonePara
{
	float lCF, rCF;
	int nChan;
	int sf;
};

class gammaToneFilter
{		
	float rPart[4], iPart[4];
		
	float dR, dI, ddR, ddI;

	float gain, twoPiT, z, f1, f2;

	int delay;

public:
	float cf, bw;

	gammaToneFilter(float cf, float bc, int sf);

	~gammaToneFilter(){;}

	void initFilter(void);
	
	void oneStep(float input);

	void filtering(float *input, int sigLength, float *response);
	
	//void filtering(float *input, int sigLength, float *response, float *insf);
};

class gammaToneFilterBank
{	
	int memIndi[MAX_CHANNEL];

	float lowerCF, upperCF;

	gammaToneFilter *gf[MAX_CHANNEL];
	
	float HzToERBRate(float Hz){ return( 21.4*log10(Hz*0.00437 + 1.0) ); }
	
	float ERBRateToHz(float ERBRate){ return( (pow(10, ERBRate/21.4) - 1) / 0.00437 ); }


public:
	int numberChannel;
	
	int sf;

	float *response[MAX_CHANNEL];
	
	gammaToneFilterBank(gammaTonePara x);

	~gammaToneFilterBank();

	void filtering(float *input, int sigLength, int chan);
};