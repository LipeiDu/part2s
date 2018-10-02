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
#include "parameters.h"
#include "ICsource_kernel.cu"
#include <H5Cpp.h>
#include <H5File.h>

using namespace std;


void readInParameters(struct parameters &params)
{
    char dummyChar[255];
    int dummyInt;
    float dummyFloat;

    FILE *fileIn;
    char fname[255];
    sprintf(fname, "parameter.dat");
    fileIn = fopen(fname, "r");

    if (fileIn == NULL)
    {
        printf("Couldn't open parameter.dat . Using default values!\n");
    }
    else
    {
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NEV = dummyInt;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.SIGMA = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.SIGMAN = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DELTA_TAU = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.TAUFORM = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.T0 = dummyFloat;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NX = dummyInt;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NY = dummyInt;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NN = dummyInt;
        fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
        params.NT = dummyInt;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DT = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DX = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DY = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.DN = dummyFloat;
        fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
        params.NORM = dummyFloat;
    }

    fclose(fileIn);
}


int main()
{

  ////////////////////////////////////////////////////////////////////////////
  //                            Initialize parameters                       //
  ////////////////////////////////////////////////////////////////////////////


  //declare parameters struct
  struct parameters params;

  // default values
  params.NEV = 1;
  params.SIGMA = 1.0;
  params.SIGMAN = 0.5;
  params.DELTA_TAU = 0.5;
  params.TAUFORM = 0.2;
  params.T0 = 0.5;
  params.NX = 131;
  params.NY = 131;
  params.NN = 61;
  params.NT = 80;
  params.DT = 0.05;
  params.DX = 0.15;
  params.DY = 0.15;
  params.DN = 0.15;
  params.NORM = 1.0;

  // read in chosen parameters from parameters.dat if such a file exists
  readInParameters(params);

  params.NTOT = (params.NX * params.NY * params.NN);

  // pass values
  float tauform = params.TAUFORM;
  int Nx = params.NX;
  int Ny = params.NY;
  int Nn = params.NN;
  int Nt = params.NT;
  int Ntot = params.NTOT;
  //float t0 = params.T0;
  //float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float dn = params.DN;
  float norm = params.NORM;

  ////////////////////////////////////////////////////////////////////////////
  //                             Converting coordinates                     //
  ////////////////////////////////////////////////////////////////////////////

  printf("SECTION of converting into Milne...\n");

  //==========================================================================
  // read in the particle list from UrQMD; get the info. in Milne.

  ifstream infile1("test.f16");

  FILE *outfile1;
  char filname[255];
  sprintf(filname, "output/Milne.dat");
  outfile1 = fopen(filname, "w");

  float r_i[4], p_i[4];

  for (int j=0; j<4; ++j)
  {
    r_i[j] = 0;
    p_i[j] = 0;
  }
  float m_i = 0;
  //float tform_i = 0;
  float b_i = 0;

  int Npart = 0; // total number of particles

  while (infile1 >> r_i[0] >> r_i[1] >> r_i[2] >> r_i[3] >> p_i[0] >> p_i[1] >> p_i[2] >> p_i[3] >> m_i)
  {

    b_i = 0.0;

    // gamma factor of particle i

    float gamma_i = p_i[0]/m_i;

    // calculate the final postion of each particle after the formation time

    float tform = tauform * gamma_i;
    float tformE = tform/p_i[0];

    r_i[0] = r_i[0] + tform;
    r_i[1] = r_i[1] + p_i[1] * tformE;
    r_i[2] = r_i[2] + p_i[2] * tformE;
    r_i[3] = r_i[3] + p_i[3] * tformE;

    // transfer into Milne coordinates

    float rm_i[4];

    rm_i[0] = sqrt(r_i[0]*r_i[0]-r_i[3]*r_i[3]); // tau
    rm_i[1] = r_i[1]; // x
    rm_i[2] = r_i[2]; // y
    rm_i[3] = 0.5 * log((r_i[0]+r_i[3])/(r_i[0]-r_i[3])); // eta_s

    float rapidity = 0.5 * log((p_i[0]+p_i[3])/(p_i[0]-p_i[3])); // rapidity
    float mT_i = sqrt(m_i*m_i+p_i[1]*p_i[1]+p_i[2]*p_i[2]); // mT

    float pm_i[4];

    pm_i[0] = mT_i * cosh(rapidity-rm_i[3]); // p_tau
    pm_i[1] = p_i[1]; // px
    pm_i[2] = p_i[2]; // py
    pm_i[3] = mT_i * sinh(rapidity-rm_i[3]); // p_eta. Caution, different definition

    // p^eta is set to 0
    pm_i[3] = 0;

    // write Milne in output file

    fprintf(outfile1,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",rm_i[0],rm_i[1],rm_i[2],rm_i[3],pm_i[0],pm_i[1],pm_i[2],pm_i[3],m_i,b_i);

    Npart++;
  }

  fclose(outfile1);

  printf("Total number of particles is %d.\n", Npart);

  //==========================================================================
  // read in the particle list in Milne

  //float ri[N][4], pi[N][4], mi[N], bi[N];
  float r0[Npart], r1[Npart], r2[Npart], r3[Npart];
  float p0[Npart], p1[Npart], p2[Npart], p3[Npart];
  float mi[Npart], bi[Npart];

  for (int i = 0; i < Npart; ++i)
  {
    r0[i] = 0.0;
    r1[i] = 0.0;
    r2[i] = 0.0;
    r3[i] = 0.0;
    p0[i] = 0.0;
    p1[i] = 0.0;
    p2[i] = 0.0;
    p3[i] = 0.0;
    mi[i] = 0.0;
    bi[i] = 0.0;
  }

  FILE *MFile;
  char fname[255];
  sprintf(fname, "output/Milne.dat");
  MFile = fopen(fname, "r");

  if (MFile == NULL) printf("The particle list in Milne couldn't be opened...\n");
  else
  {
    fseek(MFile,0L,SEEK_SET);

    for (int i = 0; i < Npart; ++i) fscanf(MFile,"%e %e %e %e %e %e %e %e %e %e", &r0[i], &r1[i], &r2[i], &r3[i], &p0[i], &p1[i], &p2[i], &p3[i], &mi[i], &bi[i]);
  }
  fclose(MFile);

  ////////////////////////////////////////////////////////////////////////////
  //                             Source terms                               //
  ////////////////////////////////////////////////////////////////////////////

  printf("Calculating source terms...\n");
  cudaDeviceSynchronize();
  cudaError_t err;

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    printf("Error at very beginning: %s\n", cudaGetErrorString(err));
    err = cudaSuccess;
  }

  //declare and allocate device arrays to hold particle info from UrQMD
  float *p0_d, *p1_d, *p2_d, *p3_d;
  float *r0_d, *r1_d, *r2_d, *r3_d;
  float *mi_d, *bi_d;
  cudaMalloc( (void**) &p0_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &p1_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &p2_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &p3_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &r0_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &r1_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &r2_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &r3_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &mi_d, Npart * sizeof(float) );
  cudaMalloc( (void**) &bi_d, Npart * sizeof(float) );

  //declare and allocate device source term arrays

  float *Sb_d, *St_d, *Sx_d, *Sy_d, *Sn_d;
  float *Ttt_d, *Ttx_d, *Tty_d, *Ttn_d, *Txx_d;
  float *Txy_d, *Txn_d, *Tyy_d, *Tyn_d, *Tnn_d;

  cudaMalloc( (void**) &Sb_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &St_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sx_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sy_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sn_d, Ntot * sizeof(float) );
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

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    printf("Error in device memory allocation: %s\n", cudaGetErrorString(err));
    err = cudaSuccess;
  }

  // copy input arrays from host to device
  cudaMemcpy( p0_d, p0, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( p1_d, p1, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( p2_d, p2, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( p3_d, p3, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( r0_d, r0, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( r1_d, r1, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( r2_d, r2, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( r3_d, r3, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( mi_d, mi, Npart * sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy( bi_d, bi, Npart * sizeof(float), cudaMemcpyHostToDevice );

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    printf("Error in a cudaMemcpy: %s\n", cudaGetErrorString(err));
    err = cudaSuccess;
  }

  //zero the device source arrays first
  cudaMemset( Sb_d, 0.0, Ntot * sizeof(float));
  cudaMemset( St_d, 0.0, Ntot * sizeof(float));
  cudaMemset( Sx_d, 0.0, Ntot * sizeof(float));
  cudaMemset( Sy_d, 0.0, Ntot * sizeof(float));
  cudaMemset( Sn_d, 0.0, Ntot * sizeof(float));

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

  //host arrays for source terms
  float *Sb, *St, *Sx, *Sy, *Sn;
  Sb = (float *)calloc( Ntot, sizeof(float) );
  St = (float *)calloc( Ntot, sizeof(float) );
  Sx = (float *)calloc( Ntot, sizeof(float) );
  Sy = (float *)calloc( Ntot, sizeof(float) );
  Sn = (float *)calloc( Ntot, sizeof(float) );

  float *Ttt, *Ttx, *Tty, *Ttn, *Txx;
  float *Txy, *Txn, *Tyy, *Tyn, *Tnn;
  Ttt = (float *)calloc( Ntot, sizeof(float) );
  Ttx = (float *)calloc( Ntot, sizeof(float) );
  Tty = (float *)calloc( Ntot, sizeof(float) );
  Ttn = (float *)calloc( Ntot, sizeof(float) );
  Txx = (float *)calloc( Ntot, sizeof(float) );
  Txy = (float *)calloc( Ntot, sizeof(float) );
  Txn = (float *)calloc( Ntot, sizeof(float) );
  Tyy = (float *)calloc( Ntot, sizeof(float) );
  Tyn = (float *)calloc( Ntot, sizeof(float) );
  Tnn = (float *)calloc( Ntot, sizeof(float) );

  //an array to hold all info for all the source terms compressed to 1d for hdf5 writer
  //float *Sall;
  //Sall = (float *)calloc( 5*Ntot, sizeof(float) );

  FILE *sourcefile, *tmunufile;
  //char source_fname[255];
  char finame[255], filename[255];

  //loop over time steps, calling kernel for each and writing to file
  for (int n = 1; n < Nt+1; ++n)
  {
    if ( (n % 10) == 1 ) printf("Calculating source term for n = %d of %d\n", n, Nt);
    int it = n-1;
    //sprintf(source_fname, "%s%d.h5", "output/Sources", n);
    //printf("Launching kernel...\n");

    source_kernel<<< grids, threads >>>(Npart, it,
                          p0_d, p1_d, p2_d, p3_d,
                          r0_d, r1_d, r2_d, r3_d,
                          mi_d, bi_d,
                          Sb_d, St_d, Sx_d, Sy_d, Sn_d,
                          Ttt_d, Ttx_d, Tty_d, Ttn_d,
                          Txx_d, Txy_d, Txn_d, Tyy_d,
                          Tyn_d, Tnn_d, params);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      printf("Error in source kernel: %s\n", cudaGetErrorString(err));
      err = cudaSuccess;
    }
    //else printf("Kernel finished, copying back to host...\n");
    //now copy results from device to host
    cudaMemcpy( Sb, Sb_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( St, St_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Sx, Sx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Sy, Sy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Sn, Sn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );

    cudaMemcpy( Ttt, Ttt_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Ttx, Ttx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Tty, Tty_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Ttn, Ttn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Txx, Txx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Txy, Txy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Txn, Txn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Tyy, Tyy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Tyn, Tyn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
    cudaMemcpy( Tnn, Tnn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );

    //compress all data into the 1d array to pass to hdf5 writer
    /*for (int is = 0; is < Ntot; is++)
    {
      Sall[is] = norm * Sb[is];
      Sall[Ntot + is] = norm * St[is];
      Sall[2 * Ntot + is] = norm * Sx[is];
      Sall[3 * Ntot + is] = norm * Sy[is];
      Sall[4 * Ntot + is] = norm * Sn[is];
    }*/
    
    //printf("Writing source terms to file...\n\n");
    //H5::H5File file(source_fname, H5F_ACC_TRUNC);
    // dataset dimensions
    //hsize_t dimsf[4];
    //dimsf[0] = Nx;
    //dimsf[1] = Ny;
    //dimsf[2] = Nn;
    //dimsf[3] = 5;

    //H5::DataSpace dataspace(4, dimsf);
    //H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    //H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    //dataset.write(Sall, H5::PredType::NATIVE_FLOAT);

    //FOR TESTING write ascii files

    sprintf(finame, "%s%d.dat", "output/Sources", n);
    sourcefile = fopen(finame, "w");
    sprintf(filename, "%s%d.dat", "output/Tmunu", n);
    tmunufile = fopen(filename, "w");

    for (int i = 0; i < Nx; ++i)
    {
      for (int j = 0; j < Ny; ++j)
      {
        for (int k = 0; k < Nn; ++k)
        {
          //float tau = t0 + ((float)n - 1.0) * dt;
          float x = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
          float y = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
          float eta = ((float)k - ((float)Nn - 1.0)/2.0) * dn;

          int s = i + j * (Nx) + k * (Nx * Ny);
          fprintf(sourcefile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, norm * St[s], norm * Sx[s], norm * Sy[s], norm * Sn[s], norm * Sb[s]);
          fprintf(tmunufile, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", x, y, eta, norm * Ttt[s], norm * Ttx[s], norm * Tty[s], norm * Ttn[s], norm * Txx[s], norm * Txy[s], norm * Txn[s], norm * Tyy[s], norm * Tyn[s], norm * Tnn[s]);
        } // for (int k )
      } //for (int j)
    } //for (int i )

    fclose(sourcefile);
    fclose(tmunufile);

    //FOR TESTING write ascii files

  } // for (int n )

  printf("Freeing memory\n");
  //clean up
  free(Sb);
  free(St);
  free(Sx);
  free(Sy);
  free(Sn);

  free(Ttt);
  free(Ttx);
  free(Tty);
  free(Ttn);
  free(Txx);
  free(Txy);
  free(Txn);
  free(Tyy);
  free(Tyn);
  free(Tnn);

  cudaFree(p0_d);
  cudaFree(p1_d);
  cudaFree(p2_d);
  cudaFree(p3_d);
  cudaFree(r0_d);
  cudaFree(r1_d);
  cudaFree(r2_d);
  cudaFree(r3_d);
  cudaFree(mi_d);
  cudaFree(bi_d);
  cudaFree(Sb_d);
  cudaFree(St_d);
  cudaFree(Sx_d);
  cudaFree(Sy_d);
  cudaFree(Sn_d);

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

  printf("Done. Goodbye! \n");
}
