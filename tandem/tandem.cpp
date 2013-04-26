#include "voicedMask.h"

// globals
double *Input;
gammaToneFilterBank *AudiPery;
voicedMask *TGroup;

// functions
void initialData(int nChan=128, double sf=16000, double theta_p=.5);

int computeData(char *filename);

void SaveOutput(FILE *fp, int frame);

void initialData(int nChan, double sf, double theta_p)
{
	featurePara x;
	gammaTonePara y;

	y.lCF = 50;	 y.rCF = 8000;
	y.nChan = nChan; y.sf = sf;
	
	x.gtP = y;
	x.bP1 = 50; x.bP2 = 450; x.bPTs = 20;
	x.theta_p = theta_p;
	
	AudiPery = new gammaToneFilterBank(x.gtP);
	
	TGroup = new voicedMask(x);
}

double* readSignalAsciiFloat(char *filename, int *numSamples1)
{
	FILE *inFile=0;
	double *signal, *signal_out;
	int numSamples=0;
	double Ftmp=0;
	int i;

	if ((inFile=fopen(filename,"r"))==NULL) {
		printf("Can not open input file %s.\n", filename);
		exit(-1);
	}

	signal=(double*)malloc(sizeof(double)*MAX_SIG_LENGTH);

	while (fscanf(inFile,"%lf",&Ftmp)!=EOF ) {
		signal[numSamples] = (double)(Ftmp);
		numSamples+=1;
	}

	fclose(inFile);

	fprintf(stderr,"processing %d samples...\n",numSamples);
	*numSamples1=numSamples;

	double sumE=0;
	for(int n=0; n<numSamples; n++)
		sumE += signal[n]*signal[n];
	sumE = double(sqrt(sumE/double(numSamples)));

	signal_out=(double*)malloc(sizeof(double)*(numSamples+1));
	for (i=0;i<numSamples;i++)
	{
		signal_out[i] = signal[i];
		signal_out[i] *= 1000/sumE;	// normalize to 60 dB average SPL
	}
	free(signal);
	return signal_out;
}

int computeData(char *filename)
{
	int sigLength, chan;

	Input=readSignalAsciiFloat(filename, &sigLength);
		
	printf("%d samples, %d frames\n", sigLength, sigLength/(TGroup->fs/100));

	TGroup->numFrame=sigLength/(TGroup->fs/100);
	
	TGroup->newPitchPara();

	for(int chan=0; chan<TGroup->numberChannel; chan++)
		AudiPery->filtering(Input, sigLength, chan);
	
	TGroup->computeFeature(AudiPery, sigLength, 1);

	TGroup->dtmPitchMask();

	delete [] Input;

	return(sigLength/(TGroup->fs/100));
}

void SaveOutput(char *outfn)
{
	char filename[255];
	FILE *fp1, *fp2;

	sprintf(filename, "%s.%d.pitch.dat", outfn, TGroup->numberChannel);
	fp1=fopen(filename, "w");

	sprintf(filename, "%s.%d.mask.dat", outfn, TGroup->numberChannel);
	fp2=fopen(filename, "w");

	fprintf(fp1, "%d %d\n", TGroup->numContour, TGroup->numFrame);

	for(int n=0; n<TGroup->numContour; n++)
	{
		for(int frame=0; frame<TGroup->numFrame; frame++)
		{
			fprintf(fp1, "%d ", TGroup->Pitch[n].value[frame]);
			
			if(TGroup->Pitch[n].value[frame]>0)
			{
				for(int chan=0; chan<TGroup->numberChannel; chan++)
					fprintf(fp2, "%6.4f ", TGroup->Pitch[n].mProb[frame].value[chan]);
				fprintf(fp2, "\n");
			}
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);
	fclose(fp2);
}

void SaveCross(char *fn_cross, char *fn_evCross)
{
	char filename1[255], filename2[255];
        FILE *fp1, *fp2;

	sprintf(filename1, "%s.%d", fn_cross, TGroup->numberChannel);
	sprintf(filename2, "%s.%d", fn_evCross, TGroup->numberChannel);
        fp1=fopen(filename1, "w");
        fp2=fopen(filename2, "w");

        for(int frame=0; frame<TGroup->numFrame; frame++)
        {
                for(int chan=0; chan<TGroup->numberChannel; chan++)
		{
			fprintf(fp1, "%6.4f ", TGroup->corrLgm[frame].cross[chan]);
			fprintf(fp2, "%6.4f ", TGroup->corrLgm[frame].crossEv[chan]);
                }
                fprintf(fp1, "\n");
		fprintf(fp2, "\n");
        }

        fclose(fp1);
        fclose(fp2);
}

void SaveEng(char *fn_eng)
{
	char filename[255];
        FILE *fp;

	sprintf(filename, "%s.%d", fn_eng, TGroup->numberChannel);
        fp=fopen(filename, "w");

        for(int frame=0; frame<TGroup->numFrame; frame++)
        {
                for(int chan=0; chan<TGroup->numberChannel; chan++)
                {
                        fprintf(fp, "%6.4f ", TGroup->corrLgm[frame].acf[chan][0]);
                }
                fprintf(fp, "\n");
        }
        fclose(fp);
}

int main(int argc, char *argv[])
{	
	if (argc!=6 && argc!=8 && argc!=9)
	{
		printf("Usage: tandem nChan sf theta_p in out \n");
		printf("       tandem nChan sf theta_p in out cross evCross \n");
		printf("       tandem nChan sf theta_p in out cross evCross eng \n");
		return 1;
	}

	char *filename;
	char *outfn;

	// initialization
	filename = argv[4];
	outfn = argv[5];
	initialData(atoi(argv[1]), atof(argv[2]), atof(argv[3]));
			
	// main program
	int numFrame=computeData(filename);

	// output	
	SaveOutput(outfn);

	if (argc==8 || argc==9)
		SaveCross(argv[6], argv[7]);

	if (argc==9)
		SaveEng(argv[8]);
		
	
	delete AudiPery;
	delete TGroup;

	return 0;
}
