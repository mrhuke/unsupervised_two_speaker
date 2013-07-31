#include "segment.h"

float *Input;

segment *TSeg;

mScaleInten *scaleInten;

gammaToneFilterBank *AudiPery;

segScale sScale;

void initialData(void);

int computeData(char *filename, char *filename2, char *filename3);

int readInput(char *filename);

void deterMineFile(FILE *fp, int count);

void SaveOutput(FILE *fp, int frame);

void initialData(int nChan)
{
	scalePara x;
	gammaTonePara y;

	y.lCF = 50;	y.rCF = 8000;
	y.nChan = nChan;	y.sf = 20000;

	x.gtP = y;

	x.lP1 = 30;	x.lP2 = 50;

	x.scale = 6;

	x.timeScale1[0] = 4; x.timeScale2[0] = 14;
	x.timeScale1[1] = 6; x.timeScale2[1] = 16;
	x.timeScale1[2] = 8; x.timeScale2[2] = 18;
	x.timeScale1[3] = 10; x.timeScale2[3] = 20;
	x.timeScale1[4] = 12; x.timeScale2[4] = 22;
	x.timeScale1[5] = 14; x.timeScale2[5] = 24;

	x.downRate = 50;

	TSeg = new segment(x.gtP.nChan);

	scaleInten = new mScaleInten(x);

	AudiPery = new gammaToneFilterBank(y);

	sScale.scale=4;

	// different scales for different filterbank
    switch (nChan)
	{
		case 64:
			sScale.timeScale[0]=0; sScale.freqScale[0]=3;
			sScale.timeScale[1]=2; sScale.freqScale[1]=3;
			sScale.timeScale[2]=5; sScale.freqScale[2]=2;
			sScale.timeScale[3]=5; sScale.freqScale[3]=2;
			break;
		case 128:
			sScale.timeScale[0]=0; sScale.freqScale[0]=6;
			sScale.timeScale[1]=2; sScale.freqScale[1]=6;
			sScale.timeScale[2]=5; sScale.freqScale[2]=6;
			sScale.timeScale[3]=5; sScale.freqScale[3]=2;
			break;
	}

	//sScale.timeScale[4]=4; sScale.freqScale[4]=6;
	//sScale.timeScale[5]=5; sScale.freqScale[5]=6;
	//sScale.timeScale[6]=5; sScale.freqScale[6]=4;
	//sScale.timeScale[7]=5; sScale.freqScale[7]=2;
	//sScale.timeScale[8]=5; sScale.freqScale[8]=1;
}

int readInput(char *filename)
{
	FILE *fp;
	fp = fopen(filename, "rb");

	fseek(fp, 0, SEEK_END);
	int sigLength = ftell(fp)/sizeof(float);

	Input = new float[sigLength];

	fseek(fp, 0, SEEK_SET);
	fread(Input, sizeof(float), sigLength, fp);
	fclose(fp);

	float sumE=0;
	for(int n=0; n<sigLength; n++)
		sumE += Input[n]*Input[n];

	sumE = float(sqrt(sumE/float(sigLength)));

	for(int n=0; n<sigLength; n++)
		Input[n] *= 1000/sumE;

	return(sigLength);
}

int readInput2(char *filename)
{
	FILE *fp;
	float f;
	int sigLength;

	fp = fopen(filename, "r");

	sigLength=0;
	while(!feof(fp))
	{
		fscanf(fp, "%f\n", &f);
		sigLength++;
	}
	fclose(fp);

	Input = new float[sigLength];

	fp = fopen(filename, "r");
	sigLength=0;
	while(!feof(fp))
	{
		fscanf(fp, "%f\n", &Input[sigLength]);
		sigLength++;
	}
	fclose(fp);

	float sumE=0;
	for(int n=0; n<sigLength; n++)
		sumE += Input[n]*Input[n];

	sumE = float(sqrt(sumE/float(sigLength)));

	for(int n=0; n<sigLength; n++)
		Input[n] *= 1000/sumE;

	return(sigLength);
}

void readCross(char *filename2, char *filename3, int numFrame)
{
	int chan, frame;
	FILE *fp;

	fp=fopen(filename2, "rb");
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<TSeg->numberChannel; chan++){
			if (chan==TSeg->numberChannel-1)
				fscanf(fp, "%f\n", &TSeg->cross[chan][frame]);
			else
				fscanf(fp, "%f ", &TSeg->cross[chan][frame]);
		}

	fclose(fp);

	fp=fopen(filename3, "rb");
	for(frame=0; frame<numFrame; frame++)
		for(chan=0; chan<TSeg->numberChannel; chan++){
			if (chan==TSeg->numberChannel-1)
				fscanf(fp, "%f\n", &TSeg->crossEv[chan][frame]);
			else
				fscanf(fp, "%f ", &TSeg->crossEv[chan][frame]);
		}

	fclose(fp);
}

int computeData(char *filename, char *filename2, char *filename3)
{
	int sigLength = readInput2(filename);
	int numFrame = sigLength/(AudiPery->sf/100);

	char drn[255];

	scalePara x;
	x.scale = 6;

	x.timeScale1[0] = 4; x.timeScale2[0] = 14;
	x.timeScale1[1] = 6; x.timeScale2[1] = 16;
	x.timeScale1[2] = 8; x.timeScale2[2] = 18;
	x.timeScale1[3] = 10; x.timeScale2[3] = 19;
	x.timeScale1[4] = 12; x.timeScale2[4] = 22;
	x.timeScale1[5] = 14; x.timeScale2[5] = 24;

	TSeg->nSample = scaleInten->nSample = sigLength/(scaleInten->downRate);

	for(int chan=0; chan<AudiPery->numberChannel; chan++)
	{
		AudiPery->filtering(Input, sigLength, chan);
	}

	scaleInten->smooth(AudiPery, sigLength, 1);

	TSeg->newSegPara();

	readCross(filename2, filename3, numFrame);

	TSeg->ofFront(scaleInten, sScale);

	TSeg->frontToSeg();

	delete [] Input;

	return(sigLength/(AudiPery->sf/100));
}

void SaveOutput(char *filename)
{
	FILE *fp;
	fp=fopen(filename, "wb");

	for(int frame=0; frame<(TSeg->nSample/TSeg->factor); frame++)
		for(int chan=0; chan<TSeg->numberChannel; chan++)
				fwrite(&TSeg->seg[chan][frame], sizeof(int), 1, fp);

	fclose(fp);
}

int main(int argc, char *argv[])
{
	if (argc!=6){
			printf("Usage: segmentation nChan in.dat in.corr in.evCorr out\n");
			exit(0);
	}

	int count=0;
	char infile1[255], infile2[255], infile3[255], outfile[255];

	sprintf(infile1, argv[2]);
	sprintf(infile2, argv[3]);
	sprintf(infile3, argv[4]);

	initialData(atoi(argv[1]));

	int numFrame=computeData(infile1, infile2, infile3);

	sprintf(outfile, argv[5]);

	SaveOutput(outfile);

	delete TSeg;
	delete scaleInten;
	delete AudiPery;
}