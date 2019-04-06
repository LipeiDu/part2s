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
#include "source_kernel.cu"
#include <H5Cpp.h>
#include <H5File.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////
//                            Parameters reader                           //
////////////////////////////////////////////////////////////////////////////

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
        params.GAMMA_MAX = dummyFloat;
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
    }

    fclose(fileIn);
}


int main()
{
  printf("\n\n\n");
  printf("++++++++++++++Program started+++++++++++++++\n");

  ////////////////////////////////////////////////////////////////////////////
  //                            Initialize parameters                       //
  ////////////////////////////////////////////////////////////////////////////


  // declare parameters struct
  struct parameters params;

  // default values
  params.NEV = 1;
  params.SIGMA = 1.0;
  params.GAMMA_MAX = 2.0;
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
  float t0 = params.T0;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float dn = params.DN;

  ////////////////////////////////////////////////////////////////////////////
  //                            Process all event files                     //
  ////////////////////////////////////////////////////////////////////////////
    
  // In this part, loop over all events; for each event, calculate the coordinates after formation time,
  // calculate the center and angle, recenter and rotate; write the coordinates in Milen in AllSets.dat;
    
  int Npart = 0; // total number of particles in all event files to be averaged over
  int Nbtot = 0;
    
  FILE *allsetfile;
  char allsetname[255];
  sprintf(allsetname, "%s.dat", "output/AllSets"); // All set.dat will be written in this file after recentering and rotation
  allsetfile = fopen(allsetname, "w");

  // loop over all events
  for (int iev=1; iev<params.NEV+1; ++iev)
  {
      
      printf("**********Processing Set%d.dat**********\n", iev);

      FILE *eventfile;
      char eventname[255];
      sprintf(eventname, "%s%d.dat", "Set",iev);
      eventfile = fopen(eventname, "r");

      ////////////////////////////////////////////////////////////////////////////
      //                             Converting coordinates                     //
      ////////////////////////////////////////////////////////////////////////////
      
      // read in the particle list from UrQMD; get the info. in Milne.
      
      printf("Converting into Milne...\n");
      
      FILE *milnefile;
      char milnename[255];
      sprintf(milnename, "%s%d.dat", "output/Milne", iev);
      milnefile = fopen(milnename, "w");
      
      float r_i[4], p_i[4];
      float rm_i[4];
      
      for (int j=0; j<4; ++j)
      {
          r_i[j] = 0;
          p_i[j] = 0;
          rm_i[j] = 0;
      }
      
      float m_i = 0;
      float tform_i = 0;
      float b_i = 0;
      float g_i = 0;// gamma factor
      
      int NpartEv = 0; // total number of particles in this event
      
      // to get the center of the energy distribution
      float Etotx = 0.0;
      float Etoty = 0.0;
      float Etot = 0.0;

      if(eventfile==NULL){
        printf("Set%d.dat could not be opened...\n", iev);
      }
      else
      {
        fseek(eventfile,0L,SEEK_SET);
        fscanf(eventfile,"%*[^\n]%*c");//Skip the header line

        while(!feof(eventfile))
        {

            // read in particle list from UrQMD
            fscanf(eventfile, "%f %f %f %f %f %f %f %f %f %f %f\n", &r_i[0], &r_i[1], &r_i[2], &r_i[3], &p_i[0], &p_i[1], &p_i[2], &p_i[3], &m_i, &tform_i, &b_i);

            // gamma factor of particle i
            g_i = p_i[0]/m_i;

            // calculate the final postion of each particle after the formation time

            float tform = tauform * g_i;
            float tformE = tform/p_i[0];

            r_i[0] = r_i[0] + tform;
            r_i[1] = r_i[1] + p_i[1] * tformE;
            r_i[2] = r_i[2] + p_i[2] * tformE;
            r_i[3] = r_i[3] + p_i[3] * tformE;

            // transfer into Milne coordinates

            rm_i[0] = sqrt(r_i[0]*r_i[0]-r_i[3]*r_i[3]); // tau
            rm_i[1] = r_i[1]; // x
            rm_i[2] = r_i[2]; // y
            rm_i[3] = 0.5 * log((r_i[0]+r_i[3])/(r_i[0]-r_i[3]+1.e-30)); // eta_s


            // write Milne in output file (WARNING: nothing is in Milne in this file, 04/05/2019)
            // skip particles outside light cone
            if ( isnan(rm_i[0]) || isnan(rm_i[3]) )
            {
              printf("*Warning* : found particle outside light cone (excluding it...)\n");
            }
            else
            {
              fprintf(milnefile,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
                      r_i[0],r_i[1],r_i[2],r_i[3],p_i[0],p_i[1],p_i[2],p_i[3],m_i,g_i,b_i);

              NpartEv++;
              Nbtot = Nbtot + b_i;

              // center of the energy distribution
              Etot = Etot + p_i[0];
              Etotx =  Etotx + p_i[0] * r_i[1];
              Etoty =  Etoty + p_i[0] * r_i[2];
            }
        }

        fclose(eventfile);
        fclose(milnefile);

        printf("Total number of particles in Set%d.dat is %d.\n", iev, NpartEv);
      }
      
      Npart = Npart + NpartEv;

      
      ////////////////////////////////////////////////////////////////////////////
      //                             Find CM and plane angle                    //
      ////////////////////////////////////////////////////////////////////////////
      
      // The center has almost been found at this point, but the angle needs to be calculated with respect to the center.
      // So the particle list needs to be looped over again

      // Section a: center
      
      float CMx = 0.0;
      float CMy = 0.0;

      CMx = Etotx/Etot;
      CMy = Etoty/Etot;

      printf("Center of Mass is: CMx= %.3f, CMy= %.3f.\n", CMx, CMy);

      for (int j=0; j<4; ++j)
      {
        rm_i[j] = 0;
        p_i[j] = 0;
      }
      
      m_i = 0;
      g_i = 0;
      b_i = 0;

      // to get the center of the energy distribution
      float psi = 0.0;
      float avgxy = 0.0;
      float avgy2x2 = 0.0;
      
      // Section b: participant plane angle
      
      milnefile = fopen(milnename, "r");
      
      if(milnefile==NULL){
          printf("Milne%d.dat could not be opened...\n", iev);
      }
      else{
          fseek(milnefile,0L,SEEK_SET);
          
          while(!feof(milnefile))
          {
              // write particle info, x in Milne, p in Cartesian
              fscanf(milnefile, "%f %f %f %f %f %f %f %f %f %f %f\n", &rm_i[0], &rm_i[1], &rm_i[2], &rm_i[3], &p_i[0], &p_i[1], &p_i[2], &p_i[3], &m_i, &g_i, &b_i);
              
              avgxy =  avgxy + p_i[0] * (rm_i[1]-CMx) * (rm_i[2]-CMy);
              avgy2x2 =  avgy2x2 + p_i[0] * ((rm_i[1]-CMx)*(rm_i[1]-CMx) - (rm_i[2]-CMy)*(rm_i[2]-CMy));
          }
      }
      
      // participant plane angle
      psi = 0.5 * atan (2 * avgxy/(avgy2x2+1.e-30));
      printf("Participant plane angle is: %.3f rad or %.3f degree.\n", psi, psi*180/3.1415);
      
      fclose(milnefile);


      ////////////////////////////////////////////////////////////////////////////
      //                             Recenter and rotate                        //
      ////////////////////////////////////////////////////////////////////////////

      float r0 = 0.0;
      float r1 = 0.0;
      float r2 = 0.0;
      float r3 = 0.0;
      float p0 = 0.0;
      float p1 = 0.0;
      float p2 = 0.0;
      float p3 = 0.0;
      float mi = 0.0;
      float gi = 0.0;
      float bi = 0.0;
    
      milnefile = fopen(milnename, "r");
      
      if(milnefile==NULL){
          printf("Milne%d.dat could not be opened...\n", iev);
      }
      else{
          fseek(milnefile,0L,SEEK_SET);
          
          while(1)
          {
              fscanf(milnefile,"%e %e %e %e %e %e %e %e %e %e %e", &r0, &r1, &r2, &r3, &p0, &p1, &p2, &p3, &mi, &gi, &bi);
              if(feof(milnefile)) break;
              
              // recenter and rotate
              
              r1 = r1 - CMx;
              r2 = r2 - CMy;
              float xp =  r1 * cos(psi) + r2 * sin(psi);
              float yp = -r1 * sin(psi) + r2 * cos(psi);
              
              float p1p =  p1 * cos(psi) + p2 * sin(psi);
              float p2p = -p1 * sin(psi) + p2 * cos(psi);
              
              fprintf(allsetfile,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",r0,xp,yp,r3,p0,p1p,p2p,p3,mi,gi,bi);
          }
      }
      
      fclose(milnefile);
  }
    
  fclose(allsetfile);
    
  printf("***************************************\n");
  printf("Total number of particles in All Sets is %d, net baryon number is %d.\n", Npart, Nbtot);

  ////////////////////////////////////////////////////////////////////////////
  //                             All particles Milne                        //
  ////////////////////////////////////////////////////////////////////////////

  // read in all the particles in all events in Cartesian

  float r0[Npart], r1[Npart], r2[Npart], r3[Npart];
  float p0[Npart], p1[Npart], p2[Npart], p3[Npart];
  float mi[Npart], gi[Npart], bi[Npart];

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
    gi[i] = 0.0;
    bi[i] = 0.0;
  }

  allsetfile = fopen(allsetname, "r");

  if (allsetfile == NULL)
    printf("The particle list in AllSets.dat couldn't be opened...\n");
  else
  {
    fseek(allsetfile,0L,SEEK_SET);
    for (int i = 0; i < Npart; ++i)
    {
        fscanf(allsetfile,"%e %e %e %e %e %e %e %e %e %e %e", &r0[i], &r1[i], &r2[i], &r3[i], &p0[i], &p1[i], &p2[i], &p3[i], &mi[i], &gi[i], &bi[i]);
    }
  }
  fclose(allsetfile);

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
  float *mi_d, *gi_d, *bi_d;
    
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

  //declare and allocate device source term arrays
  float *Sb_d, *St_d, *Sx_d, *Sy_d, *Sn_d;
  cudaMalloc( (void**) &Sb_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &St_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sx_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sy_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sn_d, Ntot * sizeof(float) );


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
  cudaMemcpy( gi_d, gi, Npart * sizeof(float), cudaMemcpyHostToDevice );
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

  //an array to hold all info for all the source terms compressed to 1d for hdf5 writer
  float *Sall;
  Sall = (float *)calloc( 5*Ntot, sizeof(float) );
    
  ////////////////////////////////////////////////////////////////////////////
  //                             Calculation all time steps                 //
  ////////////////////////////////////////////////////////////////////////////

  FILE *sourcefile;
  char source_fname[255];
  char finame[255];
    
  float Norm = 0.0;// test the normalization of the kernel, add up all time steps

  //loop over time steps, calling kernel for each and writing to file
  for (int n = 1; n < Nt+1; ++n)
  {
    if ( (n % 10) == 1 ) printf("Calculating source term for n = %d of %d\n", n, Nt);
      
    int it = n-1;
    float tau = t0 + ((float)it) * dt;
    
    printf("tau = %f\n", tau);
      
    sprintf(source_fname, "%s%d.h5", "output/Sources", n);

    source_kernel<<< grids, threads >>>(Npart, it, p0_d, p1_d, p2_d, p3_d, r0_d, r1_d, r2_d, r3_d, mi_d, gi_d, bi_d, Sb_d, St_d, Sx_d, Sy_d, Sn_d, params);

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

    //compress all data into the 1d array to pass to hdf5 writer
    for (int is = 0; is < Ntot; is++)
    {
      Sall[is] = Sb[is];
      Sall[Ntot + is] = St[is];
      Sall[2 * Ntot + is] = Sx[is];
      Sall[3 * Ntot + is] = Sy[is];
      Sall[4 * Ntot + is] = Sn[is];
    }

    //printf("Writing source terms to file...\n\n");
    H5::H5File file(source_fname, H5F_ACC_TRUNC);
      
    // dataset dimensions
    hsize_t dimsf[4];
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    dimsf[2] = Nn;
    dimsf[3] = 5;

    H5::DataSpace dataspace(4, dimsf);
    H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    dataset.write(Sall, H5::PredType::NATIVE_FLOAT);

    //FOR TESTING write ascii files
    
    sprintf(finame, "%s%d.dat", "output/Sources", n);
    sourcefile = fopen(finame, "w");
      
    for (int i = 0; i < Nx; ++i)
    {
      for (int j = 0; j < Ny; ++j)
      {
        for (int k = 0; k < Nn; ++k)
        {
          
          float x   = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
          float y   = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
          float eta = ((float)k - ((float)Nn - 1.0)/2.0) * dn;

          int s = i + j * (Nx) + k * (Nx * Ny);

          fprintf(sourcefile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, St[s], Sx[s], Sy[s], Sn[s], Sb[s]);
            
            Norm = Norm + Sb[s] * tau;
        } // for (int k )
      } //for (int j)
    } //for (int i )
    fclose(sourcefile);
      
    printf("Numerical Nb = %lf\n",Norm * dt * dx * dy * dn);
    //Norm = 0.0;
    //FOR TESTING write ascii files
  } // for (int n )
    
  Norm = Norm * dt * dx * dy * dn;
  printf("Numerical Nb = %lf\n",Norm);
    
  ////////////////////////////////////////////////////////////////////////////
  //                             Clean up                                   //
  ////////////////////////////////////////////////////////////////////////////
    
  printf("Freeing memory\n");
  
  free(Sb);
  free(St);
  free(Sx);
  free(Sy);
  free(Sn);

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

  printf("Done. Goodbye! \n");
}
