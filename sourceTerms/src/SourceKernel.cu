
#include "ParameterReader.cpp"

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
    
    float dxsqd = d0 * d0 - d1 * d1 - d2 * d2 - d3 * d3;
    float udotx = u0 * d0 - u1 * d1 - u2 * d2 - u3 * d3;
    float dist4d = dxsqd - udotx * udotx;
    
    if (dist4d <= 9 * sigma2 && udotx <= 4 * delta_tau) // if the particle's contribution is important
    {
        float exponent = dist4d * SigInv;
        float exponentiation = exp(exponent);
        
        float uxnorm = udotx * d_tauInv;
        float ch = cosh(uxnorm);
        float chInv = 1/ch;
        
        kernel = gratio * chInv * chInv * exponentiation; // [1], only cosh^2 * exp * gamma term

        if (isnan(kernel)){
            printf("Kernel is nan for this particle, set to 0...\n");
            kernel = 0.0;
        }
        
        *kernelT = gratio * exponentiation;
        
        if (isnan(*kernelT)){
            printf("KernelT is nan for this particle, set to 0...\n");
            *kernelT = 0.0;
        }

     } // contribution check
    
     return kernel;
}


__global__ void source_kernel(int Npart, int it, float *p0_d, float *p1_d, float *p2_d, float *p3_d, float *r0_d, float *r1_d, float *r2_d, float *r3_d, float *mi_d, float *gi_d, float *bi_d, float *Sb_d, float *St_d, float *Sx_d, float *Sy_d, float *Sn_d, float *Ttt_d, float *Ttx_d, float *Tty_d, float *Ttn_d, float *Txx_d, float *Txy_d, float *Txn_d, float *Tyy_d, float *Tyn_d, float *Tnn_d, parameters params)
{

  long int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
  long int cid = blockId * (blockDim.x * blockDim.y * blockDim.z) + (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;

  // parameters
  float sigma = params.SIGMA;
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
  
  // common factors
  float sigma2 = sigma*sigma;
  float SigInv = 1.0/(2.0*sigma2);
  float prefac = 1.0/pow(2.0*PI*sigma2,1.5); //[1/fm^3]
    
  float d_tauInv = 1.0/delta_tau;
  float prefactor = 0.5 * d_tauInv * prefac;
    
  float nevInv = 1.0/(float)nev;

  float facB = prefactor * nevInv; //[1/fm^4]
  float facT = prefactor * nevInv / hbarc; //[1/(GeV*fm^5)], momentum from UrQMD is in [GeV], but in Hydro code we need sources in [fm]
  float facTensor = prefac * nevInv / hbarc; //[1/(GeV*fm^4)]
    
  // gamma regulation
  float gmax2 = gmax * gmax;
  float vmax = sqrt(1-1/gmax2);
    
  // proper time
  float tau = t0 + ((float)it) * dt;
  float tauInv = 1/tau;

  
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
            
            float pm0 = cheta * p0_d[m] - sheta * p3_d[m];
            float pm3 = (-sheta * p0_d[m] + cheta * p3_d[m]) * tauInv;

            Sb = Sb + kernel * bi_d[m]; // [1]
            St = St + kernel * pm0;     // [GeV]
            Sx = Sx + kernel * p1_d[m]; // [GeV]
            Sy = Sy + kernel * p2_d[m]; // [GeV]
            Sn = Sn + kernel * pm3;     // [GeV/fm]
            
#ifdef INITIAL_TENSOR
            float ptauInv = 1 / p0_d[m]; // p0_d[m] is pt, not p^tau
            Ttt = Ttt + kernelT * ptauInv * pm0 * pm0; // kernelT is unitless
            Ttx = Ttx + kernelT * ptauInv * pm0 * p1_d[m]; // [GeV]
            Tty = Tty + kernelT * ptauInv * pm0 * p2_d[m];
            Ttn = Ttn + kernelT * ptauInv * pm0 * pm3;
            Txx = Txx + kernelT * ptauInv * p1_d[m] * p1_d[m];
            Txy = Txy + kernelT * ptauInv * p1_d[m] * p2_d[m];
            Txn = Txn + kernelT * ptauInv * p1_d[m] * pm3;
            Tyy = Tyy + kernelT * ptauInv * p2_d[m] * p2_d[m];
            Tyn = Tyn + kernelT * ptauInv * p2_d[m] * pm3;
            Tnn = Tnn + kernelT * ptauInv * pm3 * pm3; // [GeV/fm^2] = [1/GeV] * [GeV/fm] * [GeV/fm]
#endif
        } //for (int m = 0; m < N; ++m)

        //==========================================================================
        // Write the source terms to arrays
      
        Sb_d[cid] = facB * Sb; // [1/fm^4] = [1/fm^4] * [1]
        St_d[cid] = facT * St; // [1/fm^5] = [1/(fm^5*GeV)] * [GeV]
        Sx_d[cid] = facT * Sx; // [1/fm^5] = [1/(fm^5*GeV)] * [GeV]
        Sy_d[cid] = facT * Sy; // [1/fm^5] = [1/(fm^5*GeV)] * [GeV]
        Sn_d[cid] = facT * Sn; // [1/fm^6] = [1/(fm^5*GeV)] * [GeV/m]
#ifdef INITIAL_TENSOR
        Ttt_d[cid] = facTensor * Ttt; // facTensor = [1/(GeV*fm^4)]
        Ttx_d[cid] = facTensor * Ttx; // [1/fm^4] = [1/(GeV*fm^4)] * [GeV]
        Tty_d[cid] = facTensor * Tty;
        Ttn_d[cid] = facTensor * Ttn;
        Txx_d[cid] = facTensor * Txx;
        Txy_d[cid] = facTensor * Txy;
        Txn_d[cid] = facTensor * Txn;
        Tyy_d[cid] = facTensor * Tyy;
        Tyn_d[cid] = facTensor * Tyn;
        Tnn_d[cid] = facTensor * Tnn; // [1/fm^6] = [1/(GeV*fm^4)] * [GeV/fm^2]
#endif
   } //if (cid < Ntot)
}
