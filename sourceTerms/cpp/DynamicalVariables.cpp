#include "ParameterReader.cpp"

////////////////////////////////////////////////////////////////////////////
// arrays
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

float *Sb, *St, *Sx, *Sy, *Sn;
float *Sall;

float **stressTensor, **shearTensor, **flowVelocity;
float *energyDensity, *pressure, *temperature, *bulkPressure;
float **baryonCurrent, **baryonDiffusion, *baryonDensity;

float *stressAll;
float *shearAll;
float *primaryAll;
float *baryonAll;


void allocateMemory(int Ntot){
    
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

void freeMemory(){
    free(Sb);
    free(St);
    free(Sx);
    free(Sy);
    free(Sn);
    free(Sall);

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
    
    //baryon section
    for(int i = 0; i < 4; i++)
        free(baryonCurrent[i]);
    free(baryonCurrent);
    
    for(int i = 0; i < 4; i++)
        free(baryonDiffusion[i]);
    free(baryonDiffusion);
    
    free(baryonDensity);
    
    // for hdf5
    free(stressAll);
    free(shearAll);
    free(primaryAll);
    free(baryonAll);
#endif
}
