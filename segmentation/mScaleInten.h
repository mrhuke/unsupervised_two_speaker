#include "filter.h"

struct scalePara
{
	gammaTonePara gtP;

	int downRate;

	int scale;

	float timeScale1[MAX_SCALE];
	float timeScale2[MAX_SCALE];

	float lP1, lP2;
};

class mScaleInten
{
	filter *lowPass;
	filter *scalePass[MAX_SCALE];

	int memIndi;

	void newAmpScale(void);
	
	void deleteAmpScale(void);


public:
	int nSample, numberChannel;

	int fs, downRate;

	int scale;

	float *ampScale[MAX_CHANNEL*MAX_SCALE];
	
	mScaleInten(scalePara x);
	
	~mScaleInten(void);

	void smooth(gammaToneFilterBank *AudiPery, int sigLength, int echo);

	void readInten(int i, int j, int snr, scalePara, int sigLength, char *drn);
	void writeInten(int i, int j, int snr, scalePara, int sigLength, char *drn);
};
