
#include <math.h>
#include "ParameterReader.cpp"
//#define MILNE_KERNEL

__host__ __device__ float KernelCartesian(float r0, float r1, float r2, float r3, float ri0, float ri1, float ri2, float ri3, float pi0, float pi1, float pi2, float pi3, float mi, float gi, float gmax, float vmax, float sigma2, float delta_tau, float d_tauInv, float SigInv, float *kernelT){
    
    // ********************************************************
    // if the contraction factor is too large, rescale it
    // ********************************************************
    
    float mInv = 1/mi;
    float u1 = pi1 * mInv;
    float u2 = pi2 * mInv;
    float u3 = pi3 * mInv;
    float u0 = sqrt(1 + u1*u1 + u2*u2 + u3*u3);
    
    float gratio = gmax / gi;
    
    if (gi > gmax)
    {
        float gi2 = gi * gi;
        float vi = sqrt(1-1/gi2);
        float vratio = vmax / vi;
        float ratio = gratio * vratio;
        
        u0 = u0 * gratio;
        u1 = u1 * ratio;
        u2 = u2 * ratio;
        u3 = u3 * ratio;
        
        gi = u0;
    }
    
    // ********************************************************
    // calculate the kernel
    // ********************************************************
    
    float kernel = 0.0;
    
    //convert (tau, eta) into (t, z)
    float t = r0 * cosh(r3);
    float z = r0 * sinh(r3);
    
    float d0 = t  - ri0;
    float d1 = r1 - ri1;
    float d2 = r2 - ri2;
    float d3 = z  - ri3;
    
    float distTem = d0 / gi;
    float dr2 = d1 * d1 + d2 * d2 + d3 * d3;
    float udr = u1 * d1 + u2 * d2 + u3 * d3;
    float distSpa2 = dr2 + udr * udr;
    
    if (distSpa2 <= 9 * sigma2 && distTem <= 3 * delta_tau) // if the particle's contribution is important
    {
        float expSpace = - distSpa2 * SigInv;
        float expS = exp(expSpace);
        
        float distTem2 = distTem * distTem;
        float expTem = exp(- distTem2 * d_tauInv * d_tauInv / 2);
        
        // kernel for dynamical sources after taking derivatives
        kernel = expS * expTem; // [1], not including the factors

        if (isnan(kernel)){
            printf("Kernel is nan for this particle, set to 0...\n");
            kernel = 0.0;
        }
        
        // kernel for tensors before taking derivatives; fluid tensor, not particle tensor
        float erfExp = 0.5 * ( erf(distTem * d_tauInv / 1.41421) + 1 );
        
        *kernelT = erfExp * expS;
        
        if (isnan(*kernelT)){
            printf("KernelT is nan for this particle, set to 0...\n");
            *kernelT = 0.0;
        }

     } // contribution check
    
     return kernel;
}

__host__ __device__ float KernelMilne(float r1, float r2, float r3, float rm1, float rm2, float rm3, float sigma, float sigman, float SigInv, float SignInv){

    float kernel = 0.0;
    
    float distn = fabs(r3 - rm3);
    float distn2 = distn * distn;
    
    if (distn < 4 * sigman)//if it's not far away in eta direction
    {
        float ddx = fabs(r1 - rm1);
        float ddy = fabs(r2 - rm2);
        
        float disttr2 = ddx*ddx + ddy*ddy;
        float disttr = sqrt(disttr2);
        
        if (disttr < 4 * sigma)//if the particle is not far away in the transverse plane
        {
            float dist = -(disttr2 * SigInv + distn2 * SignInv);
            kernel = exp(dist);
        }
    }
    
    if (isnan(kernel)){
        printf("Kernel is nan for this particle, set to 0...\n");
        kernel = 0.0;
    }
    
    return kernel;
}


__global__ void source_kernel(int Npart, int it, float *p0_d, float *p1_d, float *p2_d, float *p3_d, float *r0_d, float *r1_d, float *r2_d, float *r3_d, float *mi_d, float *gi_d, float *bi_d, float *Sb_d, float *St_d, float *Sx_d, float *Sy_d, float *Sn_d, float *Ttt_d, float *Ttx_d, float *Tty_d, float *Ttn_d, float *Txx_d, float *Txy_d, float *Txn_d, float *Tyy_d, float *Tyn_d, float *Tnn_d, parameters params)
{

    long int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
    long int cid = blockId * (blockDim.x * blockDim.y * blockDim.z) + (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;

    // parameters
    float sigma = params.SIGMA;
    float sigman = params.SIGMAN;
    float delta_tau = params.DELTA_TAU;
    float gmax = params.GAMMA_MAX;
    float t0 = params.T0;
    int nev = params.NEV;
    int Ntot = params.NTOT;
    int Nx = params.NX;
    int Ny = params.NY;
    int Nn = params.NN;
    float dt = params.DT;
    float dx = params.DX;
    float dy = params.DY;
    float dn = params.DN;
    
    float nevInv = 1.0/(float)nev;
  
    // common factors
    float sigma2 = sigma * sigma;
    float sigma3 = sigma2 * sigma;
    float SigInv = 1.0 / (2.0 * sigma2);
    float spaceGaussianFactor = 1.0 / pow(2.0*PI,1.5) / sigma3; //[1/fm^3]
    
    float d_tauInv = 1.0/delta_tau;
    float timeGaussianFactor = 1.0 / sqrt(2.0*PI) / delta_tau; //[1/fm]
    
    float combinedFactor = spaceGaussianFactor * timeGaussianFactor; //[1/fm^4]
    
    float facSource = combinedFactor * nevInv; //[1/fm^4]
    float facTensor = spaceGaussianFactor * nevInv; //[1/(fm^3)]
    
    // gamma regulation
    float gmax2 = gmax * gmax;
    float vmax = sqrt(1-1/gmax2);
    
    // proper time
    float tau = t0 + ((float)it) * dt;
    float tauInv = 1/tau;
    
    // smearing in Milne
    float prefacMilne = 1/(pow(2.0*PI,1.5) * sigma2 * sigman) * tauInv;
    float SignInv =  1.0/(2.0*sigman*sigman);
    float facTensorMilne = prefacMilne * nevInv;
  
    //==========================================================================
    // loop over all cells (x, y, eta)
    
    if (cid < Ntot)
    {
      
        // reconstruct indices manually using
        int k = cid / (Nx * Ny);
        int j = ( cid - (k * Nx * Ny) ) / Nx;
        int i = cid - (k * Nx * Ny) - (j * Nx);

        // space-time where to calculate sources
        float r0 = tau;
        float r1 = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
        float r2 = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
        float r3 = ((float)k - ((float)Nn - 1.0)/2.0) * dn;
      
        // coordinate transformation matrix element
        float cheta = cosh(r3); // cosh(eta_s)
        float sheta = sinh(r3); // sinh(eta_s)
      
        // calculate the source terms for rid ijk and write them in the output file
        float Sb = 0.0;
        float St = 0.0;
        float Sx = 0.0;
        float Sy = 0.0;
        float Sn = 0.0;
      
#ifdef INITIAL_TENSOR
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
#endif

        // ****************************************************************************
        // OPTION I: Smearing in Cartesian with Lorentz contraction
        // ****************************************************************************
      
#ifndef MILNE_KERNEL
        //==========================================================================
        // loop over all particles

        for (int m = 0; m < Npart; ++m)
        {

            // ********************************************************
            // calculate the kernel
            // ********************************************************
            
            float kernel = 0.0;
            float kernelT = 0.0;
            
            kernel = KernelCartesian(r0, r1, r2, r3, r0_d[m], r1_d[m], r2_d[m], r3_d[m], p0_d[m], p1_d[m], p2_d[m], p3_d[m],
                                     mi_d[m], gi_d[m], gmax, vmax, sigma2, delta_tau, d_tauInv, SigInv, &kernelT);
            
            // ********************************************************
            // calculate the source
            // ********************************************************
            
            float pm0 = cheta * p0_d[m] - sheta * p3_d[m]; // [1]
            float pm3 = (-sheta * p0_d[m] + cheta * p3_d[m]) * tauInv; // [1/fm]
            
            float u0Inv = 1 / gi_d[m]; // 1/u0

            Sb = Sb + kernel * u0Inv * bi_d[m]; // [1]
            St = St + kernel * u0Inv * pm0;     // [1/fm]
            Sx = Sx + kernel * u0Inv * p1_d[m]; // [1/fm]
            Sy = Sy + kernel * u0Inv * p2_d[m]; // [1/fm]
            Sn = Sn + kernel * u0Inv * pm3;     // [1/fm^2]
            
#ifdef INITIAL_TENSOR
            float ptauInv = 1 / p0_d[m]; // p0_d[m] is pt, not p^tau
            Ttt = Ttt + kernelT * ptauInv * pm0 * pm0; // [1/fm]
            Ttx = Ttx + kernelT * ptauInv * pm0 * p1_d[m];
            Tty = Tty + kernelT * ptauInv * pm0 * p2_d[m];
            Ttn = Ttn + kernelT * ptauInv * pm0 * pm3;
            Txx = Txx + kernelT * ptauInv * p1_d[m] * p1_d[m];
            Txy = Txy + kernelT * ptauInv * p1_d[m] * p2_d[m];
            Txn = Txn + kernelT * ptauInv * p1_d[m] * pm3;
            Tyy = Tyy + kernelT * ptauInv * p2_d[m] * p2_d[m];
            Tyn = Tyn + kernelT * ptauInv * p2_d[m] * pm3;
            Tnn = Tnn + kernelT * ptauInv * pm3 * pm3;
#endif
        } //for (int m = 0; m < N; ++m)

        //==========================================================================
        // Write the source terms to arrays
      
        Sb_d[cid] = facSource * Sb; // [1/fm^4] = [1/fm^4] * [1]
        St_d[cid] = facSource * St; // [1/fm^5] = [1/fm^4] * [1/fm]
        Sx_d[cid] = facSource * Sx; // [1/fm^5] = [1/fm^4] * [1/fm]
        Sy_d[cid] = facSource * Sy; // [1/fm^5] = [1/fm^4] * [1/fm]
        Sn_d[cid] = facSource * Sn; // [1/fm^6] = [1/fm^4] * [1/fm^2]
#ifdef INITIAL_TENSOR
        Ttt_d[cid] = facTensor * Ttt; // facTensor = [1/fm^3]
        Ttx_d[cid] = facTensor * Ttx; // [1/fm^4] = [1/fm^3] * [1/fm]
        Tty_d[cid] = facTensor * Tty;
        Ttn_d[cid] = facTensor * Ttn;
        Txx_d[cid] = facTensor * Txx;
        Txy_d[cid] = facTensor * Txy;
        Txn_d[cid] = facTensor * Txn;
        Tyy_d[cid] = facTensor * Tyy;
        Tyn_d[cid] = facTensor * Tyn;
        Tnn_d[cid] = facTensor * Tnn;
#endif
      
      // ****************************************************************************
      // OPTION II: Smearing in Milne
      // ****************************************************************************
#else
      for (int m = 0; m < Npart; ++m)
      {
          
          float kernelT = 0.0;
          
          float rm0 = sqrt(r0_d[m] * r0_d[m] - r3_d[m] * r3_d[m]);
          float rm3 = 0.5 * log((r0_d[m] + r3_d[m]) / (r0_d[m] - r3_d[m] + 1.e-30));
          float pm0 = cheta * p0_d[m] - sheta * p3_d[m];
          float pm3 = (-sheta * p0_d[m] + cheta * p3_d[m]) * tauInv;
          
          kernelT = KernelMilne(r1, r2, r3, r1_d[m], r2_d[m], rm3, sigma, sigman, SigInv, SignInv);

          float ptauInv = 1 / pm0; // pm0 is p^tau
          
          Ttt = Ttt + kernelT * ptauInv * pm0 * pm0; // kernelT is unitless
          Ttx = Ttx + kernelT * ptauInv * pm0 * p1_d[m];
          Tty = Tty + kernelT * ptauInv * pm0 * p2_d[m];
          Ttn = Ttn + kernelT * ptauInv * pm0 * pm3;
          Txx = Txx + kernelT * ptauInv * p1_d[m] * p1_d[m];
          Txy = Txy + kernelT * ptauInv * p1_d[m] * p2_d[m];
          Txn = Txn + kernelT * ptauInv * p1_d[m] * pm3;
          Tyy = Tyy + kernelT * ptauInv * p2_d[m] * p2_d[m];
          Tyn = Tyn + kernelT * ptauInv * p2_d[m] * pm3;
          Tnn = Tnn + kernelT * ptauInv * pm3 * pm3; //
      } //for (int m = 0; m < N; ++m)
      
      Sb_d[cid] = 0.0;
      St_d[cid] = 0.0;
      Sx_d[cid] = 0.0;
      Sy_d[cid] = 0.0;
      Sn_d[cid] = 0.0;

      Ttt_d[cid] = facTensorMilne * Ttt;
      Ttx_d[cid] = facTensorMilne * Ttx;
      Tty_d[cid] = facTensorMilne * Tty;
      Ttn_d[cid] = facTensorMilne * Ttn;
      Txx_d[cid] = facTensorMilne * Txx;
      Txy_d[cid] = facTensorMilne * Txy;
      Txn_d[cid] = facTensorMilne * Txn;
      Tyy_d[cid] = facTensorMilne * Tyy;
      Tyn_d[cid] = facTensorMilne * Tyn;
      Tnn_d[cid] = facTensorMilne * Tnn;
#endif
   } //if (cid < Ntot)
}
