////////////////////////////////////////////////////////////////////////////
//                                  Parameters                            //
////////////////////////////////////////////////////////////////////////////

#define hbarc 0.1973
#define PI 3.14159

// UrQMD events
#define nev 1 // # of UrQMD events

// smearing width
#define sigma 1.0 // width in transverse plane
#define sigman 0.5 // width in eta direction

// thermalization time and width
#define delta_tau 0.5
#define tauform 0.2

// source terms related
#define t0 0.5 // time of the 1st souce term
#define Nx 131
#define Ny 131
#define Nn 61
#define Nt 80
#define Ntot (Nx * Ny * Nn)
#define dt 0.05
#define dx 0.15
#define dy 0.15
#define dn 0.15
