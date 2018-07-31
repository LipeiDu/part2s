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

using namespace std;

int main()
{

  float SigInv = 1.0/(2.0*sigma*sigma);
  float SignInv = 1.0/(2.0*sigman*sigman);
  float d_tauInv = 1.0/delta_tau;
  float nevInv = 1.0/(float)nev;
  float hbarcNevInv = nevInv/hbarc;

  float prefac = 1.0/(2.0 * (2.0*PI*sigma*sigma) * sqrt(2.0*PI*sigman*sigman));
  float prefactor = d_tauInv * prefac;
  float facN = prefactor * nevInv;
  float facHN = prefactor * hbarcNevInv;

  ////////////////////////////////////////////////////////////////////////////
  //                             Converting coordinates                     //
  ////////////////////////////////////////////////////////////////////////////

  printf("SECTION of converting into Milne...\n");

  //==========================================================================
  // read in the particle list from UrQMD; get the info. in Milne.

  ifstream infile1("Set.dat");

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
  float tform_i = 0;
  float b_i = 0;

  int Npart = 0; // total number of particles

  while (infile1 >> r_i[0] >> r_i[1] >> r_i[2] >> r_i[3] >> p_i[0] >> p_i[1] >> p_i[2] >> p_i[3] >> m_i >> tform_i >> b_i)
  {
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

  if(MFile==NULL){
    printf("The particle list in Milne couldn't be opened...\n");
  }
  else{
    fseek(MFile,0L,SEEK_SET);
    for(int i = 0; i < Npart; ++i)
    fscanf(MFile,"%e %e %e %e %e %e %e %e %e %e", &r0[i], &r1[i], &r2[i], &r3[i], &p0[i], &p1[i], &p2[i], &p3[i], &mi[i], &bi[i]);
  }

  fclose(MFile);

  ////////////////////////////////////////////////////////////////////////////
  //                             Source terms                               //
  ////////////////////////////////////////////////////////////////////////////

  printf("Calculating source terms...\n");

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
  cudaMalloc( (void**) &Sb_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &St_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sx_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sy_d, Ntot * sizeof(float) );
  cudaMalloc( (void**) &Sn_d, Ntot * sizeof(float) );

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

  //launch the cuda kernel that computes the source terms
  // need enough threads to cover Ntot total number of spacetime points
  int threadsPerBlock = 256;
  int numBlocks = (Ntot / threadsPerBlock) + 1;
  source_kernel<<<numBlocks, threadsPerBlock>>>(Npart,
                        p0_d, p1_d, p2_d, p3_d,
                        r0_d, r1_d, r2_d, r3_d,
                        mi_d, bi_d,
                        Sb_d, St_d, Sx_d, Sy_d, Sn_d);

  //now copy results from device to host
  float Sb[Ntot], St[Ntot], Sx[Ntot], Sy[Ntot], Sn[Ntot];
  cudaMemcpy( Sb, Sb_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
  cudaMemcpy( St, St_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
  cudaMemcpy( Sx, Sx_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
  cudaMemcpy( Sy, Sy_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );
  cudaMemcpy( Sn, Sn_d, Ntot * sizeof(float), cudaMemcpyDeviceToHost );

  //write results to file

  for (int n = 1; n < Nt+1; ++n)
  {
    int it = n-1;

    FILE *sourcefile;
    char finame[255];
    sprintf(finame, "%s%d.dat", "output/Sources", n);
    sourcefile = fopen(finame, "w");

    for (int i = 0; i < Nx; ++i)
    {
      for (int j = 0; j < Ny; ++j)
      {
        for (int k = 0; k < Nn; ++k)
        {
          float tau = t0 + ((float)n - 1.0) * dt;
          float x = ((float)i - ((float)Nx + 1.0)/2.0) * dx;
          float y = ((float)j - ((float)Ny + 1.0)/2.0) * dy;
          float eta = ((float)k - ((float)Nn + 1.0)/2.0) * dn;

          int s = it * (Nx * Ny * Nn) + i * (Ny * Nn) + j * (Nn) + k;
          fprintf(sourcefile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, St[s], Sx[s], Sy[s], Sn[s], Sb[s]);
        } // for (int k )
      } //for (int j)
    } //for (int i )

    fclose(sourcefile);
  } // for (int n )
  //clean up
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
}
