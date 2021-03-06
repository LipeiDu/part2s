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

    if (fileIn == NULL) printf("Couldn't open parameter.dat . Using default values!\n");
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
    }
    fclose(fileIn);
}

int main()
{

  ////////////////////////////////////////////////////////////////////////////
  //                            Initialize parameters                       //
  ////////////////////////////////////////////////////////////////////////////

  // declare parameters struct
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

  float max_x = (float)(Nx-1) * dx / 2.0;
  float max_y = (float)(Ny-1) * dy / 2.0;
  float max_eta = (float)(Nn-1) * dn / 2.0;

  ////////////////////////////////////////////////////////////////////////////
  //                             Converting coordinates                     //
  ////////////////////////////////////////////////////////////////////////////

    FILE *outfile;
    char filname[255];
    sprintf(filname, "output/EventSummary.dat");
    outfile = fopen(filname, "w");
    
  // read in all event files, each event is in a single file
    
  for(int n = 1; n < params.NEV+1; ++n){
      FILE *eventfile;
      char fname[255];
      sprintf(fname, "%s%d.dat", "input/event",n);
      eventfile = fopen(fname, "r");
      
      int nParticipant = 0;
      int event = 0;
      float impactPara = 0.0;
      
      if(eventfile==NULL){
          printf("The event file %d could not be opened...\n", n);
          //exit(-1);
      }
      else
      {
          fseek(eventfile,0L,SEEK_SET);
          
          // read and print the header line
          fscanf(eventfile,"%d%f%d", &event, &impactPara, &nParticipant);
          //printf("event=%d\t impact parameter=%f\t No. of participant=%d \n", event, impactPara, nParticipant);
          
          //printf("Converting to Milne coordinates and integrating total event energy \n");
          float total_energy_in = 0.0;
          float total_energy_out = 0.0;
          
          //==========================================================================
          // read in the particle list from UrQMD; get the info. in Milne.

          float r_i[4], p_i[4];
          for (int j = 0; j < 4; ++j)
          {
              r_i[j] = 0.0;
              p_i[j] = 0.0;
          }
          
          float m_i = 0;
          float tform_i = 0;
          float b_i = 0;
          
          int Npart = 0; // total number of particles
          
          fscanf(eventfile,"%*[^\n]%*c");//Skip the header line
              
          for(int i=0;!feof(eventfile);++i)
          {
              fscanf(eventfile, "%f %f %f %f %f %f %f %f %f %f %f\n", &r_i[0], &r_i[1], &r_i[2], &r_i[3], &p_i[0], &p_i[1], &p_i[2], &p_i[3], &m_i, &tform_i, &b_i);
              
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
              
              //check if particle is inside hydro box, if so add its energy to total inside box
              //2D box
              if (Nn == 1)
              {
                  if ( (fabs(rm_i[1]) < max_x)  && (fabs(rm_i[2]) < max_y)  && !isnan(rm_i[0])) total_energy_in += p_i[0];
                  //otherwise keep track of energy outside hydro box
                  else total_energy_out += p_i[0];
              }
              //3D box
              else
              {
                  if ( (fabs(rm_i[1]) < max_x)  && (fabs(rm_i[2]) < max_y) && (fabs(rm_i[3]) < max_eta)  && !isnan(rm_i[0])) total_energy_in += p_i[0];
                  //otherwise keep track of energy outside hydro box
                  else total_energy_out += p_i[0];
              }
              Npart++;
          }
          
          //printf("Total number of particles           :  %d.\n", Npart);
          //printf("Total energy inside hydro box       :  %f.\n", total_energy_in);
          //printf("Total energy outside hydro box      :  %f.\n", total_energy_out);
          float energy_ratio = total_energy_in / (total_energy_in + total_energy_out);
          //printf("Percentage energy intside hydro box :  %f.\n", energy_ratio);

          // write information to EventSummary.dat, event ID, impact parameter, number of participants, number of produced particles, total energy inside hydro box, Total energy outside hydro box, Percentage energy intside hydro box
          fprintf(outfile,"%d\t %.8f\t %d\t %d\t %.8f\t %.8f\t %.8f\n", event, impactPara, nParticipant, Npart, total_energy_in, total_energy_out, energy_ratio);
          
          fclose(eventfile);
        }
    }
    
    fclose(outfile);
    
    printf("Done. Goodbye! \n");
}
