////////////////////////////////////////////////////////////////////////////
//                                  Parameters                            //
////////////////////////////////////////////////////////////////////////////
#pragma once

#define hbarc 0.1973
#define PI 3.14159

//#define INITIAL_TENSOR

struct parameters
{
// UrQMD events
    int NEV;// # of UrQMD events

// smearing width
    float SIGMA; // width in transverse plane
    float SIGMAN;
    float GAMMA_MAX; // cutoff of the Lorentz contraction

// thermalization time and width
    float DELTA_TAU;
    float TAUFORM;

// source terms related
    float T0; // time of the 1st souce term
    int NX;
    int NY;
    int NN;
    int NT;
    int NTOT;
    float DT;
    float DX;
    float DY;
    float DN;
};


////////////////////////////////////////////////////////////////////////////
//                            Parameters reader                           //
////////////////////////////////////////////////////////////////////////////

void readInParameters(struct parameters &params)
{
    // default values
    params.NEV = 1;
    params.SIGMA = 1.0;
    params.SIGMAN = 1.0;
    params.GAMMA_MAX = 2.0;
    params.DELTA_TAU = 0.5;
    params.TAUFORM = 0.2;
    params.T0 = 0.5;
    params.NX = 131;
    params.NY = 131;
    params.NN = 61;
    params.NT = 80;
    params.DT = 0.05;
    params.DX = 0.15;
    params.DY = 0.15;
    params.DN = 0.15;
    
    // read in parameters from parameter.dat
    char dummyChar[255];
    int dummyInt;
    float dummyFloat;
    
    FILE *fileIn;
    char fname[255];
    sprintf(fname, "parameter.dat");
    fileIn = fopen(fname, "r");
    
    if (fileIn == NULL)
    {
        printf("Couldn't open parameter.dat . Using default values!\n");
    }
    else
    {
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NEV = dummyInt;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.SIGMA = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.SIGMAN = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.GAMMA_MAX = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DELTA_TAU = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.TAUFORM = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.T0 = dummyFloat;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NX = dummyInt;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NY = dummyInt;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NN = dummyInt;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NT = dummyInt;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DT = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DX = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DY = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DN = dummyFloat;
    }
    
    fclose(fileIn);
    
    // total cubes
    params.NTOT = (params.NX * params.NY * params.NN);
}



