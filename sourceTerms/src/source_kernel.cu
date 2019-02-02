
#include "parameters.h"

__global__ void source_kernel(int Npart, int it,
                              float *p0_d, float *p1_d, float *p2_d, float *p3_d,
                              float *r0_d, float *r1_d, float *r2_d, float *r3_d,
                              float *mi_d, float *bi_d,
                              float *Sb_d, float *St_d, float *Sx_d, float *Sy_d, float *Sn_d, parameters params)
{

  /*
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int k = threadIdx.z + blockIdx.z * blockDim.z;
  long int tid = i + j * (blockDim.x * gridDim.x) + k * (blockDim.x * gridDim.x * blockDim.y * gridDim.y);
  */

  long int blockId = blockIdx.x + blockIdx.y * gridDim.x
                  + gridDim.x * gridDim.y * blockIdx.z;
  long int tid =     blockId * (blockDim.x * blockDim.y * blockDim.z)
                  + (threadIdx.z * (blockDim.x * blockDim.y))
                  + (threadIdx.y * blockDim.x) + threadIdx.x;

  // parameters
  int Ntot = params.NTOT;
  int nev = params.NEV;
  float sigma = params.SIGMA;

  float delta_tau = params.DELTA_TAU;
  float t0 = params.T0;
  int Nx = params.NX;
  int Ny = params.NY;
  int Nn = params.NN;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float dn = params.DN;
    
  float sigma2 = sigma*sigma;
  float SigInv = 1.0/(2.0*sigma2);
  float prefac = 1.0/pow(2.0*PI*sigma2,1.5);
    
  float d_tauInv = 1.0/delta_tau;
  float prefactor = 0.5 * d_tauInv * prefac;
    
  float nevInv = 1.0/(float)nev;

  float facB = prefactor * nevInv; //[1/fm^4]
  float facT = prefactor * nevInv / hbarc; //[1/GeV*fm^5]
  
  // loop over all cells
  if (tid < Ntot)
  {
    
    // reconstruct indices manually using
    int k = tid / (Nx * Ny);
    int j = ( tid - (k * Nx * Ny) ) / Nx;
    int i = tid - (k * Nx * Ny) - (j * Nx);

    //printf("tid = %d\n", tid);
    float tau = t0 + ((float)it) * dt;
    float tau2 = tau * tau;
    float tauInv = 1/tau;

    //==========================================================================
    // calculate the source terms for rid ijk and write them in the output file

    float Sb = 0.0;
    float St = 0.0;
    float Sx = 0.0;
    float Sy = 0.0;
    float Sn = 0.0;

    // space-time where to calculate sources
    float r0 = tau;
    float r1 = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
    float r2 = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
    float r3 = ((float)k - ((float)Nn - 1.0)/2.0) * dn;

    //==========================================================================
    // loop over all particles

    for (int m = 0; m < Npart; ++m)
    {
      if ( !isnan(r0_d[m]) )//if the particle is outside the light cone, skip this particle
      {
        float b_i = bi_d[m];
          
        float d0 = r0 - r0_d[m];
        float d1 = r1 - r1_d[m];
        float d2 = r2 - r2_d[m];
        float d3 = r3 - r3_d[m];
          
        float mInv = 1/mi_d[m];
        float u0 = p0_d[m] * mInv;
        float u1 = p1_d[m] * mInv;
        float u2 = p2_d[m] * mInv;
        float u3 = p3_d[m] * mInv; // u3 doesn't have 1/tau
          
        float dxsqd = d0 * d0 - d1 * d1 - d2 * d2 - tau2 * d3 * d3;
        float udotx = u0 * d0 - u1 * d1 - u2 * d2 - u3 * d3;
        float dist4d = dxsqd - udotx * udotx;
          
        if (dist4d <= 9 * sigma2 && udotx <= 3 * delta_tau) // if the particle's contribution is important
        {
          // Smearing kernel
          float exponent = dist4d * SigInv;
          float exponentiation = exp(exponent);
            
          float uxnorm = udotx * d_tauInv;
          float ch = cosh(uxnorm);
          float chInv = 1/ch;
          float kernel = chInv * chInv * exponentiation; // [1], only the cosh^2 and exp product

          if (!isnan(kernel))// if kernel is nan for some reasons, skip this particle
          {
            Sb = Sb + kernel * b_i;     // [1]
            St = St + kernel * p0_d[m]; // [GeV]
            Sx = Sx + kernel * p1_d[m]; // [GeV]
            Sy = Sy + kernel * p2_d[m]; // [GeV]
            Sn = Sn + kernel * p3_d[m]; // [GeV], pi[m][4] is [GeV] by defination above
           }
            
         } // contribution check
       } //if ( !isnan(ri[m][0]) )
     } //for (int m = 0; m < N; ++m)

    //==========================================================================
    // Write the source terms to arrays
    Sb_d[tid] = facB * Sb; // [1/fm^4] = [1/fm^4] * [1]
    St_d[tid] = facT * St; // [1/fm^5] = [1/(fm^5*GeV)] * [GeV]
    Sx_d[tid] = facT * Sx; // [1/fm^5] = [1/(fm^5*GeV)] * [GeV]
    Sy_d[tid] = facT * Sy; // [1/fm^5] = [1/(fm^5*GeV)] * [GeV]
    Sn_d[tid] = facT * Sn * tauInv; // [1/fm^6] = [1/(fm^5*GeV)] * [GeV] * [1/fm]
      
   } //if (tid < Ntot)
}
