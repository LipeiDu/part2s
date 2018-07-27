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

using namespace std;

////////////////////////////////////////////////////////////////////////////
//                                  Parameters                            //
////////////////////////////////////////////////////////////////////////////

#define hbarc 0.1973
#define PI 3.14159

// UrQMD events
#define nev 50 // # of UrQMD events

// smearing width
#define sigma 0.5 // width in transverse plane
#define sigman 0.2 // width in eta direction

// thermalization time and width
#define delta_tau 0.5
#define tauform 0.2

// source terms related
#define t0 0.5 // time of the 1st souce term
#define Nx 135
#define Ny 135
#define Nn 1
#define Nt 1
#define dt 0.05
#define dx 0.15
#define dy 0.15
#define dn 0.15

int main()
{

    float SigInv = 1/(2*sigma*sigma);
    float SignInv = 1/(2*sigman*sigman);
    float d_tauInv = 1/delta_tau;
    float nevInv = 1/(float)nev;
    float hbarcNevInv = nevInv/hbarc;

    float prefac = 1/(2 * (2*PI*sigma*sigma) * sqrt(2*PI*sigman*sigman));
    float prefactor = d_tauInv * prefac;
    float facN = prefactor * nevInv;
    float facHN = prefactor * hbarcNevInv;

////////////////////////////////////////////////////////////////////////////
//                             Converting coordinates                     //
////////////////////////////////////////////////////////////////////////////
    
    printf("SECTION of converting into Milne...\n");

    //==========================================================================
    // read in the particle list from UrQMD; get the info. in Milne.
    
    ifstream infile1("../../Set.dat");
    
    FILE *outfile1;
    char filname[255];
    sprintf(filname, "../../output/Milne.dat");
    outfile1 = fopen(filname, "w");
    
    float r_i[4], p_i[4], m_i, tform_i, b_i;
    
    for(int j=0; j<4; ++j){
        r_i[j] = 0;
        p_i[j] = 0;
    }
    m_i = 0;
    tform_i = 0;
    b_i = 0;
    
    int N=0; // total number of particles
    
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
        
        N++;
    }
    
    fclose(outfile1);

    printf("Total number of particles is %d.\n", N);

    //==========================================================================
    // read in the particle list in Milne
    
    float ri[N][4], pi[N][4], mi[N], bi[N];
    
    for(int i=0; i<N; ++i){
        for(int j=0; j<4; ++j){
            ri[i][j] = 0;
            pi[i][j] = 0;
        }
        mi[i] = 0;
        bi[i] = 0;
    }
    
    FILE *MFile;
    char fname[255];
    sprintf(fname, "../../output/Milne.dat");
    MFile = fopen(fname, "r");
    
    if(MFile==NULL){
        printf("The particle list in Milne couldn't be opened...\n");
    }
    else{
        fseek(MFile,0L,SEEK_SET);
        for(int i=0; i<N; ++i)
            fscanf(MFile,"%e %e %e %e %e %e %e %e %e %e", & ri[i][0], & ri[i][1], & ri[i][2], & ri[i][3], & pi[i][0], & pi[i][1], & pi[i][2], & pi[i][3], & mi[i], & bi[i]);
    }
    
    fclose(MFile);
    
////////////////////////////////////////////////////////////////////////////
//                             Source terms                               //
////////////////////////////////////////////////////////////////////////////
    
    printf("SECTION of calculating source terms...\n");
    
    // Calculate the source terms
    
    float time = t0;
    
    clock_t t1, t2;
    
    for(int n = 1; n < Nt+1; ++n){
        
        //==========================================================================
        // time step l
        
        printf("Calculate source terms for the time step %d\n", n);
        
        FILE *sourcefile;
        char finame[255];
        sprintf(finame, "%s%d.dat", "../../output/Sources", n);
        sourcefile = fopen(finame, "w");
        
        float t1 = clock();
             
        float tauInv = 1/time;
        float rr[4];

        
        for(int i = 0; i < Nx; ++i){
            for(int j = 0; j < Ny; ++j){
                for(int k = 0; k < Nn; ++k){
                    
                    //==========================================================================
                    // calculate the source terms for rid ijk and write them in the output file
                    
                    float Sb = 0;
                    float St = 0;
                    float Sx = 0;
                    float Sy = 0;
                    float Sn = 0;
                    
                    rr[0] = time;
                    rr[1] = ((float)i - ((float)Nx + 1.0)/2.0) * dx;
                    rr[2] = ((float)j - ((float)Ny + 1.0)/2.0) * dy;
                    rr[3] = ((float)k - ((float)Nn + 1.0)/2.0) * dn;
                    
                    //==========================================================================
                    // loop over all particles
                    
                    for(int m = 0; m < N; ++m){
                        
                        float kernel = 0;
                        
                        if(!isnan(ri[m][0])){ //if the particle is outside the light cone, skip this particle
                            
                            b_i = bi[m];
                            
                            float distt = fabs(rr[0]-ri[m][0]);
                            
                            if(distt < 4 * delta_tau)//if it's not far away in tau
                            {
                                
                                float distn = fabs(rr[3]-ri[m][3]);
                                
                                if(distn < 3 * sigman)//if it's not far away in eta direction
                                {
                                    float ddx = fabs(rr[1]-ri[m][1]);
                                    float ddy = fabs(rr[2]-ri[m][2]);
                                    float disttrs = ddx*ddx + ddy*ddy;
                                    float disttr = sqrt(disttrs);
                                    
                                    if(disttr < 3 * sigma)//if the particle is not far away in the transverse plane
                                    {
                                        // Smearing kernel
                                        
                                        float dist = -(disttrs * SigInv + distn*distn * SignInv);
                                        float numerator = exp(dist);
                                        float delta = distt * d_tauInv;
                                        float ch = cosh(delta);
                                        
                                        kernel = 1/(ch * ch) * numerator; // [1]
                                        
                                        if(!isnan(kernel)){ // if kernel is nan for some reasons, skip this particle
                                            
                                            // pi[m][4] is [GeV] by defination above
                                            
                                            Sb = Sb + kernel * b_i; // [1]
                                            St = St + kernel * pi[m][0]; // [GeV]
                                            Sx = Sx + kernel * pi[m][1]; // [GeV]
                                            Sy = Sy + kernel * pi[m][2]; // [GeV]
                                            Sn = Sn + kernel * pi[m][3] * tauInv; // [GeV/fm] caution, definition and tau here
                                            
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    //==========================================================================
                    // Write the source terms in output files
                    
                    // [1/fm^3] facN = 1/(2 * delta_tau * (2*PI*sigma*sigma) * sqrt(2*PI*sigman*sigman))/nev;
                    // [1/(fm^4*GeV)] facHN = 1/(2 * delta_tau * (2*PI*sigma*sigma) * sqrt(2*PI*sigman*sigman))/nev/hbarc;
                    
                    Sb = facN  * Sb * tauInv; // [1/fm^4] = [1/fm^3] * [1] * [1/fm]
                    St = facHN * St * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
                    Sx = facHN * Sx * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
                    Sy = facHN * Sy * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
                    Sn = facHN * Sn * tauInv; // [1/fm^6] = [1/(fm^4*GeV)] * [GeV/fm] * [1/fm]
                    
                    //if(fabs(St) < 1.e-10) St = 0.0;
                    //if(fabs(Sx) < 1.e-10) Sx = 0.0;
                    //if(fabs(Sy) < 1.e-10) Sy = 0.0;
                    //if(fabs(Sn) < 1.e-10) Sn = 0.0;
                    //if(fabs(Sb) < 1.e-10) Sb = 0.0;
                    
                    fprintf(sourcefile,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",rr[1],rr[2],rr[3],St,Sx,Sy,Sn,Sb);
                    
                }
            }
        }
        
        fclose(sourcefile);
        
        time = t0 + n * dt;
        
        float t2 = clock();
        float delta_time = (t2 - t1)/(double)(CLOCKS_PER_SEC);
        printf("(Elapsed time: %.3f ms)\n", delta_time);
    }
    
    return 0;
}