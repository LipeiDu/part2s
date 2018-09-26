
#include "parameters.h"

__global__ void source_kernel(int Npart, int it,
                              float *p0_d, float *p1_d, float *p2_d, float *p3_d,
                              float *r0_d, float *r1_d, float *r2_d, float *r3_d,
                              float *mi_d, float *bi_d,
                              float *Sb_d, float *St_d, float *Sx_d, float *Sy_d, float *Sn_d,
                              float *Ttt_d, float *Ttx_d, float *Tty_d, float *Ttn_d,
                              float *Txx_d, float *Txy_d, float *Txn_d, float *Tyy_d,
                              float *Tyn_d, float *Tnn_d, parameters params)
{

  // components of Tmunu: Ttt, Ttx, Tty, Ttn, Txx, Txy, Txn, Tyy, Tyn, Tnn

  long int blockId = blockIdx.x + blockIdx.y * gridDim.x
                  + gridDim.x * gridDim.y * blockIdx.z;
  long int tid =     blockId * (blockDim.x * blockDim.y * blockDim.z)
                  + (threadIdx.z * (blockDim.x * blockDim.y))
                  + (threadIdx.y * blockDim.x) + threadIdx.x;

  int Ntot = params.NTOT;
  int nev = params.NEV;
  float sigma = params.SIGMA;
  float sigman = params.SIGMAN;
  float delta_tau = params.DELTA_TAU;
  float t0 = params.T0;
  int Nx = params.NX;
  int Ny = params.NY;
  int Nn = params.NN;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float dn = params.DN;
  
  if (tid < Ntot)
  {
    
    //reconstruct indices manually using
    // s = i + j * (Nx) + k * (Nx * Ny);
    int k = tid / (Nx * Ny);
    int j = ( tid - (k * Nx * Ny) ) / Nx;
    int i = tid - (k * Nx * Ny) - (j * Nx);

    //printf("tid = %d\n", tid);
    float tau = t0 + ((float)it) * dt;
    float tauInv = 1.0 / tau;
    float rr[4];

    float SigInv = 1.0/(2.0*sigma*sigma);
    float SignInv = 1.0/(2.0*sigman*sigman);
    float d_tauInv = 1.0/delta_tau;
    float nevInv = 1.0/(float)nev;
    float hbarcNevInv = nevInv/hbarc;

    float prefac = 1.0/(2.0 * (2.0*PI*sigma*sigma) * sqrt(2.0*PI*sigman*sigman));
    float prefactor = d_tauInv * prefac;
    float facN = prefactor * nevInv;
    float facHN = prefactor * hbarcNevInv;
    float hc2NevInv = nevInv/(hbarc*hbarc);

    //==========================================================================
    // calculate the source terms for rid ijk and write them in the output file

    float Sb = 0.0;
    float St = 0.0;
    float Sx = 0.0;
    float Sy = 0.0;
    float Sn = 0.0;

    float Ttt = 0.0;
    float Ttx = 0.0;
    float Tty = 0.0;
    float Ttn = 0.0;
    float Txx = 0.0;
    float Txy = 0.0;
    float Txn = 0.0;
    float Tyy = 0.0;
    float Tyn = 0.0;
    float Tnn = 0.0;

    rr[0] = tau;
    rr[1] = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
    rr[2] = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
    rr[3] = ((float)k - ((float)Nn - 1.0)/2.0) * dn;

    //==========================================================================
    // loop over all particles

    for (int m = 0; m < Npart; ++m)
    {
      if ( !isnan(r0_d[m]) )
      { //if the particle is outside the light cone, skip this particle
        float b_i = bi_d[m];
        float distt = fabs(rr[0]-r0_d[m]);
        if (distt < 4 * delta_tau)//if it's not far away in tau
        {
          float distn = fabs(rr[3]-r3_d[m]);
          if (distn < 3 * sigman)//if it's not far away in eta direction
          {
            float ddx = fabs(rr[1]-r1_d[m]);
            float ddy = fabs(rr[2]-r2_d[m]);
            float disttrs = ddx*ddx + ddy*ddy;
            float disttr = sqrt(disttrs);

            if (disttr < 3 * sigma)//if the particle is not far away in the transverse plane
            {
              // Smearing kernel
              float dist = -(disttrs * SigInv + distn*distn * SignInv);
              float numerator = exp(dist);
              float delta = distt * d_tauInv;
              float ch = cosh(delta);
              float kernel = 1.0/(ch * ch) * numerator; // [1]
              
              float mptau = mi_d[m]/p0_d[m];

              // Tmunu related
              float th = tanh((rr[0]-r0_d[m])*d_tauInv);
              float kernelT = (th + 1) * numerator;

              if ( !isnan(kernel) )
              { // if kernel is nan for some reasons, skip this particle
                // pi[m][4] is [GeV] by defination above
                Sb = Sb + kernel * mptau * b_i; // [1]
                St = St + kernel * mptau * p0_d[m]; // [GeV]
                Sx = Sx + kernel * mptau * p1_d[m]; // [GeV]
                Sy = Sy + kernel * mptau * p2_d[m]; // [GeV]
                Sn = Sn + kernel * mptau * p3_d[m] * tauInv; // [GeV/fm] caution, definition and tau here
              }

              if ( !isnan(kernelT) )
              {
                // components of Tmunu: Ttt, Ttx, Tty, Ttn, Txx, Txy, Txn, Tyy, Tyn, Tnn

                float ptauInv = 1/p0_d[m];
                Ttt = Ttt + kernelT * p0_d[m];
                Ttx = Ttx + kernelT * ptauInv * p0_d[m] * p1_d[m];
                Tty = Tty + kernelT * ptauInv * p0_d[m] * p2_d[m];
                Ttn = Ttn + kernelT * ptauInv * p0_d[m] * p3_d[m];
                Txx = Txx + kernelT * ptauInv * p1_d[m] * p1_d[m];
                Txy = Txy + kernelT * ptauInv * p1_d[m] * p2_d[m];
                Txn = Txn + kernelT * ptauInv * p1_d[m] * p3_d[m];
                Tyy = Tyy + kernelT * ptauInv * p2_d[m] * p2_d[m];
                Tyn = Tyn + kernelT * ptauInv * p2_d[m] * p3_d[m];
                Tnn = Tnn + kernelT * ptauInv * p3_d[m] * p3_d[m];
              }
            } //if (disttr < 3 * sigma)
          } //if (distn < 3 * sigman)
        } //if (distt < 4 * delta_tau)
      } //if ( !isnan(ri[m][0]) )
    } //for (int m = 0; m < N; ++m)

    //==========================================================================
    // Write the source terms to arrays
    Sb_d[tid] = facN  * Sb * tauInv; // [1/fm^4] = [1/fm^3] * [1] * [1/fm]
    St_d[tid] = facHN * St * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
    Sx_d[tid] = facHN * Sx * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
    Sy_d[tid] = facHN * Sy * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
    Sn_d[tid] = facHN * Sn * tauInv; // [1/fm^6] = [1/(fm^4*GeV)] * [GeV/fm] * [1/fm]

    float hctfac = tauInv * hc2NevInv;

    Ttt_d[tid] = prefac * Ttt * hctfac;
    Ttx_d[tid] = prefac * Ttx * hctfac;
    Tty_d[tid] = prefac * Tty * hctfac;
    Ttn_d[tid] = prefac * Ttn * hctfac;
    Txx_d[tid] = prefac * Txx * hctfac;
    Txy_d[tid] = prefac * Txy * hctfac;
    Txn_d[tid] = prefac * Txn * hctfac;
    Tyy_d[tid] = prefac * Tyy * hctfac;
    Tyn_d[tid] = prefac * Tyn * hctfac;
    Tnn_d[tid] = prefac * Tnn * hctfac;
  } //if (tid < Ntot)
}
