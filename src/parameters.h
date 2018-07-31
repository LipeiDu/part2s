////////////////////////////////////////////////////////////////////////////
//                                  Parameters                            //
////////////////////////////////////////////////////////////////////////////

#define hbarc 0.1973
#define PI 3.14159

// UrQMD events
#define nev 50 // # of UrQMD events

// smearing width
#define sigma 0.5 // width in transverse plane
#define sigman 0.2 // width in eta direction

// thermalization time and width
#define delta_tau 0.5
#define tauform 0.2

// source terms related
#define t0 0.5 // time of the 1st souce term
#define Nx 31
#define Ny 31
#define Nn 31
#define Nt 20
#define Ntot (Nt * Nx * Ny * Nn)
#define dt 0.05
#define dx 0.5
#define dy 0.5
#define dn 0.5
