////////////////////////////////////////////////////////
//               DynamicalSource.cpp                  //
//                                                    //
//         Created by Lipei Du on 7/23/18.            //
////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <math.h>
#include <string>
#include <iomanip>
#include "ParameterReader.cpp"
#include "ProcessEvents.cpp"
#include "LandauMatch.cpp"
#include "SourceKernel.cu"
#include "DynamicalVariables.cu"
#include "FileWriter.cpp"

using namespace std;

int main()
{
  printf("\n");
  printf("++++++++++++++Program started+++++++++++++++\n");

  ////////////////////////////////////////////////////////////////////////////
  //                            Initialize parameters                       //
  ////////////////////////////////////////////////////////////////////////////

  // declare parameters struct and initialization
  struct parameters params;
  readInParameters(params);

  // pass values
  float tauform = params.TAUFORM;
  int Nx = params.NX;
  int Ny = params.NY;
  int Nn = params.NN;
  int Nt = params.NT;
  int Ntot = params.NTOT;
  float t0 = params.T0;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float dn = params.DN;

  ////////////////////////////////////////////////////////////////////////////
  //                            Process all event files                     //
  ////////////////////////////////////////////////////////////////////////////
    
  // In this part, loop over all events; for each event, calculate the coordinates after formation time,
  // calculate the center and angle, recenter and rotate; write the coordinates in AllSets.dat;
    
  int Npart = 0; // total number of particles in all event files to be averaged over
  
#ifndef INITIAL_TENSOR
  int Nbtot = 0;
  eventsTreatment(params.NEV, tauform, &Nbtot, &Npart);
  printf("Total number of particles in All Sets is %d, net baryon number is %d.\n", Npart, Nbtot);
#else
  testfileTreatment(params.NEV, tauform, &Npart);
  printf("Total number of particles in test.16 is %d.\n", Npart);
#endif

  ////////////////////////////////////////////////////////////////////////////
  //                             All particles                              //
  ////////////////////////////////////////////////////////////////////////////

  // read in all the particles in all events in Cartesian

  float r0[Npart], r1[Npart], r2[Npart], r3[Npart];
  float p0[Npart], p1[Npart], p2[Npart], p3[Npart];
  float mi[Npart], gi[Npart], bi[Npart];

  readInFormedParticles(Npart, r0, r1, r2, r3, p0, p1, p2, p3, mi, gi, bi);


  ////////////////////////////////////////////////////////////////////////////
  //                             Allocate memory                            //
  ////////////////////////////////////////////////////////////////////////////

  cudaDeviceSynchronize();
  cudaError_t err;

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    printf("Error at very beginning: %s\n", cudaGetErrorString(err));
    err = cudaSuccess;
  }

  //declare and allocate device arrays
  allocateDeviceMemory(Npart, Ntot, err);
    
  //copy input arrays from host to device
  copyHostToDeviceMemory(Npart, Ntot, r0, r1, r2, r3, p0, p1, p2, p3, mi, gi, bi, err);

  allocateHostMemory(Ntot);
  
  //kernel launch parameters
  int threadsX = 16;
  int threadsY = 16;
  int threadsZ = 1;
  int blocksX = (Nx+threadsX-1)/threadsX;
  int blocksY = (Ny+threadsY-1)/threadsY;
  int blocksZ = (Nn+threadsZ-1)/threadsZ;
    
  printf("CUDA kernel parameters:\n");
  printf("dim3 grids = ( %d, %d, %d )\n", blocksX, blocksY, blocksZ);
  printf("dim3 threads = ( %d, %d, %d )\n", threadsX, threadsY, threadsZ);
    
  dim3 grids( blocksX, blocksY, blocksZ );
  dim3 threads( threadsX, threadsY, threadsZ);
    
  ////////////////////////////////////////////////////////////////////////////
  //                             Calculation all time steps                 //
  ////////////////////////////////////////////////////////////////////////////
    
  printf("***************************************\n");
  printf("Calculating source terms...\n");
    
  float Norm = 0.0;// test the normalization of the kernel, add up all time steps
  float T00total = 0.0;
  float S0total = 0.0;

  //loop over time steps, calling kernel for each and writing to file
  for (int n = 1; n < Nt+1; ++n)
  {
    if ( (n % 10) == 1 ) printf("Calculating source term for n = %d of %d\n", n, Nt);
      
    int it = n-1;
    float tau = t0 + ((float)it) * dt;
    float dV = tau * dt * dx * dy * dn;
    
    printf("tau = %f\n", tau);
      
    // conservation test
    *dT00tot = 0.0;
    *dS0tot = 0.0;
    
    source_kernel<<< grids, threads >>>(Npart, it, p0_d, p1_d, p2_d, p3_d, r0_d, r1_d, r2_d, r3_d, mi_d, gi_d, bi_d,
                                        Sb_d, St_d, Sx_d, Sy_d, Sn_d, Ttt_d, Ttx_d, Tty_d, Ttn_d, Txx_d, Txy_d, Txn_d,
                                        Tyy_d, Tyn_d, Tnn_d, dT00tot_d, dS0tot_d, params);

    copyDeviceToHostMemory(Ntot, err);
      
    T00total = T00total + (*dT00tot) * dV * hbarc * params.NEV;
    S0total = S0total + (*dS0tot) * dV * hbarc * params.NEV;
    printf("Integrated T00 is %lf, integrated S0 is %lf.\n", T00total, S0total);
      
    // Landau matching
#ifdef INITIAL_TENSOR
    solveEigenSystem(stressTensor, energyDensity, flowVelocity, pressure, tau, Nx, Ny, Nn, Ntot, dx, dy);
    calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure, Ntot, tau);
    calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor, Ntot, tau);
#endif
      
    ////////////////////////////////////////////////////////////////////////////
    //                             write output file                          //
    ////////////////////////////////////////////////////////////////////////////

    // write HDF5 file
    writeSourcesHDF5(n, Nx, Ny, Nn, Ntot, Sall, Sb, St, Sx, Sy, Sn);

    //FOR TESTING write ascii files
    writeSourcesASCII(n, Nx, Ny, Nn, Ntot, tau, dt, dx, dy, dn, Sb, St, Sx, Sy, Sn, &Norm);
      
    // write tensor file
    writeTensorsASCII(n, Nx, Ny, Nn, Ntot, tau, dt, dx, dy, dn, Sb, St, Sx, Sy, Sn, &Norm);

  } // for (int n ) time steps
    
    printf("FINAL: Integrated T00 is %lf, integrated S0 is %lf.\n", T00total, S0total);
    
  ////////////////////////////////////////////////////////////////////////////
  //                             Clean up                                   //
  ////////////////////////////////////////////////////////////////////////////
    
  printf("Freeing memory\n");
  freeMemory();

  printf("Done. Goodbye! \n");
}
