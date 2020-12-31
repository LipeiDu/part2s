#include "ParameterReader.cpp"
#include <cuda.h>

////////////////////////////////////////////////////////////////////////////
// device arrays
////////////////////////////////////////////////////////////////////////////

//declare and allocate device arrays to hold particle info from UrQMD
extern float *p0_d, *p1_d, *p2_d, *p3_d;
extern float *r0_d, *r1_d, *r2_d, *r3_d;
extern float *mi_d, *gi_d, *bi_d;

//declare and allocate device source term arrays and tensors
extern float *Sb_d, *St_d, *Sx_d, *Sy_d, *Sn_d;
extern float *Ttt_d, *Ttx_d, *Tty_d, *Ttn_d, *Txx_d, *Txy_d, *Txn_d, *Tyy_d, *Tyn_d, *Tnn_d;
extern float *Nt_d, *Nx_d, *Ny_d, *Nn_d;

////////////////////////////////////////////////////////////////////////////
// host arrays
////////////////////////////////////////////////////////////////////////////

//host arrays for source terms
extern float *Sb, *St, *Sx, *Sy, *Sn;

//an array to hold all info for all the source terms compressed to 1d for hdf5 writer
extern float *Sall;

//host arrays for tensor
extern float **stressTensor, **shearTensor, **flowVelocity;
extern float *energyDensity, *pressure, *temperature, *bulkPressure;
extern float **baryonCurrent, **baryonDiffusion, *baryonDensity;

//arrays for hdf5 writer
extern float *stressAll;
extern float *shearAll;
extern float *primaryAll;
extern float *baryonAll;

////////////////////////////////////////////////////////////////////////////
// allocation
////////////////////////////////////////////////////////////////////////////

float *p0_d, *p1_d, *p2_d, *p3_d;
float *r0_d, *r1_d, *r2_d, *r3_d;
float *mi_d, *gi_d, *bi_d;

float *Sb_d, *St_d, *Sx_d, *Sy_d, *Sn_d;
float *Ttt_d, *Ttx_d, *Tty_d, *Ttn_d, *Txx_d, *Txy_d, *Txn_d, *Tyy_d, *Tyn_d, *Tnn_d;
float *Nt_d, *Nx_d, *Ny_d, *Nn_d;

float *Sb, *St, *Sx, *Sy, *Sn;
float *Sall;

float **stressTensor, **shearTensor, **flowVelocity;
float *energyDensity, *pressure, *temperature, *bulkPressure;
float **baryonCurrent, **baryonDiffusion, *baryonDensity;

float *stressAll;
float *shearAll;
float *primaryAll;
float *baryonAll;

////////////////////////////////////////////////////////////////////////////

void allocateDeviceMemory(int Npart, int Ntot, cudaError_t err){
    
    cudaMalloc( (void**) &p0_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &p1_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &p2_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &p3_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &r0_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &r1_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &r2_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &r3_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &mi_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &gi_d, Npart * sizeof(float) );
    cudaMalloc( (void**) &bi_d, Npart * sizeof(float) );
    
    cudaMalloc( (void**) &Sb_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &St_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Sx_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Sy_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Sn_d, Ntot * sizeof(float) );
    
    //zero the device source arrays first
    cudaMemset( Sb_d, 0.0, Ntot * sizeof(float));
    cudaMemset( St_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Sx_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Sy_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Sn_d, 0.0, Ntot * sizeof(float));
    
#ifdef INITIAL_TENSOR
    //declare and allocate device tensor arrays
    cudaMalloc( (void**) &Ttt_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Ttx_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Tty_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Ttn_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Txx_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Txy_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Txn_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Tyy_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Tyn_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Tnn_d, Ntot * sizeof(float) );
    
    cudaMalloc( (void**) &Nt_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Nx_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Ny_d, Ntot * sizeof(float) );
    cudaMalloc( (void**) &Nn_d, Ntot * sizeof(float) );
    
    //zero the device tensor arrays first
    cudaMemset( Ttt_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Ttx_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Tty_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Ttn_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Txx_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Txy_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Txn_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Tyy_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Tyn_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Tnn_d, 0.0, Ntot * sizeof(float));
    
    cudaMemset( Nt_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Nx_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Ny_d, 0.0, Ntot * sizeof(float));
    cudaMemset( Nn_d, 0.0, Ntot * sizeof(float));
#endif
    
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in device memory allocation: %s\n", cudaGetErrorString(err));
        err = cudaSuccess;
    }
}

void allocateHostMemory(int Ntot){
    
    //host arrays for source terms
    Sb = (float *)calloc( Ntot, sizeof(float) );
    St = (float *)calloc( Ntot, sizeof(float) );
    Sx = (float *)calloc( Ntot, sizeof(float) );
    Sy = (float *)calloc( Ntot, sizeof(float) );
    Sn = (float *)calloc( Ntot, sizeof(float) );
    
    //an array to hold all info for all the source terms compressed to 1d for hdf5 writer
    Sall = (float *)calloc( 5*Ntot, sizeof(float) );
    
#ifdef INITIAL_TENSOR
    //host arrays for tensor
    stressTensor = (float **)calloc( 10, sizeof(float*));
    for(int i = 0; i < 10; i++)
        stressTensor[i] = (float *)calloc( Ntot, sizeof(float) );
    
    shearTensor = (float **)calloc( 10, sizeof(float*));
    for(int i = 0; i < 10; i++)
        shearTensor[i] = (float *)calloc( Ntot, sizeof(float) );
    
    flowVelocity = (float **)calloc( 4, sizeof(float*));
    for(int i = 0; i < 4; i++)
        flowVelocity[i] = (float *)calloc( Ntot, sizeof(float) );
    
    energyDensity = (float *)calloc( Ntot, sizeof(float) );
    pressure  = (float *)calloc( Ntot, sizeof(float) );
    temperature = (float *)calloc( Ntot, sizeof(float) );
    bulkPressure = (float *)calloc( Ntot, sizeof(float) );
    
    //baryon sector
    baryonCurrent = (float **)calloc( 4, sizeof(float*));
    for(int i = 0; i < 4; i++)
        baryonCurrent[i] = (float *)calloc( Ntot, sizeof(float) );
    
    baryonDiffusion = (float **)calloc( 4, sizeof(float*));
    for(int i = 0; i < 4; i++)
        baryonDiffusion[i] = (float *)calloc( Ntot, sizeof(float) );
    
    baryonDensity = (float *)calloc( Ntot, sizeof(float) );
    
    //for hdf5
    stressAll = (float *)calloc( 10*Ntot, sizeof(float) ); //Tmunu
    shearAll = (float *)calloc( 10*Ntot, sizeof(float) ); //shear stress
    primaryAll = (float *)calloc( 8*Ntot, sizeof(float) ); //energy, pressure, temperature, bulk, flow velocity
    baryonAll = (float *)calloc( 9*Ntot, sizeof(float) ); //baryon density, net current, diffusion
#endif
}

void copyHostToDeviceMemory(int Npart, int Ntot, float *r0, float *r1, float *r2, float *r3, float *p0, float *p1, float *p2, float *p3, float *mi, float *gi, float *bi, cudaError_t err){
    
    cudaMemcpy( p0_d, p0, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( p1_d, p1, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( p2_d, p2, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( p3_d, p3, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( r0_d, r0, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( r1_d, r1, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( r2_d, r2, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( r3_d, r3, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( mi_d, mi, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( gi_d, gi, Npart * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( bi_d, bi, Npart * sizeof(float), cudaMemcpyHostToDevice );
    
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in a cudaMemcpy: %s\n", cudaGetErrorString(err));
        err = cudaSuccess;
    }
}

void copyDeviceToHostMemory(int Ntot, cudaError_t err){
    
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in source kernel: %s\n", cudaGetErrorString(err));
        err = cudaSuccess;
    }
    
    //now copy results from device to host
    cudaMemcpy( Sb, Sb_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( St, St_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Sx, Sx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Sy, Sy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Sn, Sn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    
#ifdef INITIAL_TENSOR
    cudaMemcpy( stressTensor[0], Ttt_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[1], Ttx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[2], Tty_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[3], Ttn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[4], Txx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[5], Txy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[6], Txn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[7], Tyy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[8], Tyn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( stressTensor[9], Tnn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    
    cudaMemcpy( baryonCurrent[0], Nt_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( baryonCurrent[1], Nx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( baryonCurrent[2], Ny_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( baryonCurrent[3], Nn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
#endif
    
}

void freeMemory(){
    free(Sb);
    free(St);
    free(Sx);
    free(Sy);
    free(Sn);
    free(Sall);
    
    cudaFree(p0_d);
    cudaFree(p1_d);
    cudaFree(p2_d);
    cudaFree(p3_d);
    cudaFree(r0_d);
    cudaFree(r1_d);
    cudaFree(r2_d);
    cudaFree(r3_d);
    cudaFree(mi_d);
    cudaFree(gi_d);
    cudaFree(bi_d);
    cudaFree(Sb_d);
    cudaFree(St_d);
    cudaFree(Sx_d);
    cudaFree(Sy_d);
    cudaFree(Sn_d);

#ifdef INITIAL_TENSOR
    for(int i = 0; i < 10; i++)
        free(stressTensor[i]);
    free(stressTensor);
    
    for(int i = 0; i < 10; i++)
        free(shearTensor[i]);
    free(shearTensor);
    
    for(int i = 0; i < 4; i++)
        free(flowVelocity[i]);
    free(flowVelocity);
    
    free(energyDensity);
    free(pressure);
    free(temperature);
    free(bulkPressure);
    
    cudaFree(Ttt_d);
    cudaFree(Ttx_d);
    cudaFree(Tty_d);
    cudaFree(Ttn_d);
    cudaFree(Txx_d);
    cudaFree(Txy_d);
    cudaFree(Txn_d);
    cudaFree(Tyy_d);
    cudaFree(Tyn_d);
    cudaFree(Tnn_d);
    
    //baryon section
    for(int i = 0; i < 4; i++)
        free(baryonCurrent[i]);
    free(baryonCurrent);
    
    for(int i = 0; i < 4; i++)
        free(baryonDiffusion[i]);
    free(baryonDiffusion);
    
    free(baryonDensity);

    cudaFree(Nt_d);
    cudaFree(Nx_d);
    cudaFree(Ny_d);
    cudaFree(Nn_d);
    
    // for hdf5
    free(stressAll);
    free(shearAll);
    free(primaryAll);
    free(baryonAll);
#endif
}
