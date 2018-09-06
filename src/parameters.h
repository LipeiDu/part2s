////////////////////////////////////////////////////////////////////////////
//                                  Parameters                            //
////////////////////////////////////////////////////////////////////////////
#pragma once

#define hbarc 0.1973
#define PI 3.14159

struct parameters
{
// UrQMD events
    int NEV;// # of UrQMD events

// smearing width
    float SIGMA; // width in transverse plane
    float SIGMAN; // width in eta direction

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
