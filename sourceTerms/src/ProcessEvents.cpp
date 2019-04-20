
#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

////////////////////////////////////////////////////////////////////////////
//                            Process all events                          //
////////////////////////////////////////////////////////////////////////////

void eventsTreatment(int nev, double tauform, int *Nbtot, int *Npart){
    
    printf("Read UrQMD output files Set.dat with baryon number...\n");
    
    float energyTotal = 0.0;
    float massTotal = 0.0;
    
    FILE *allsetfile;
    char allsetname[255];
    sprintf(allsetname, "%s.dat", "output/AllSets"); // All set.dat will be written in this file after recentering and rotation
    allsetfile = fopen(allsetname, "w");
    
    // loop over all events
    for (int iev=1; iev<nev+1; ++iev)
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
                
                // total energy and mass
                energyTotal = energyTotal + p_i[0];
                massTotal = massTotal + m_i;
                
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
                    *Nbtot = *Nbtot + b_i;
                    
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
        
        *Npart = *Npart + NpartEv;
        
        
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
        
        if(isnan(psi)){
            printf("Participant plane angle is nan!");
            exit(-1);
        }else
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
    
    printf("Total energy of all particle is %lf, total mass is %lf.\n",energyTotal,massTotal);
    
}

////////////////////////////////////////////////////////////////////////////
//                            UrQMD test.16 file                          //
////////////////////////////////////////////////////////////////////////////

void testfileTreatment(int nev, double tauform, int *Npart){
    
    printf("Read original UrQMD output file test.16...\n");
    
    // total energy and mass
    float energyTotal = 0.0;
    float massTotal = 0.0;
    
    //==========================================================================
    // read in the particle list from UrQMD; get the info. in Milne.
    
    ifstream infile1("test.f16");
    
    FILE *outfile1;
    char filname[255];
    sprintf(filname, "output/AllSets.dat");
    outfile1 = fopen(filname, "w");
    
    float r_i[4], p_i[4];
    
    for (int j=0; j<4; ++j)
    {
        r_i[j] = 0;
        p_i[j] = 0;
    }
    
    float m_i = 0;
    float b_i = 0;
    
    while (infile1 >> r_i[0] >> r_i[1] >> r_i[2] >> r_i[3] >> p_i[0] >> p_i[1] >> p_i[2] >> p_i[3] >> m_i)
    {
        
        b_i = 0.0;
        
        // gamma factor of particle i
        
        float gamma_i = p_i[0]/m_i;
        
        // total energy and mass
        energyTotal = energyTotal + p_i[0];
        massTotal = massTotal + m_i;
        
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
        
        // write particles in output file
        
        if ( isnan(rm_i[0]) || isnan(rm_i[3]) )
        {
            printf("*Warning* : found particle outside light cone (excluding it...)\n");
        }
        else
        {
            fprintf(outfile1,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
                    r_i[0],r_i[1],r_i[2],r_i[3],p_i[0],p_i[1],p_i[2],p_i[3],m_i,gamma_i,b_i);
            
            *Npart = *Npart + 1;
        }
    }
    
    fclose(outfile1);
    
    printf("Total energy of all particle is %lf, total mass is %lf.\n",energyTotal,massTotal);

}

////////////////////////////////////////////////////////////////////////////
//                            AllSet.dat                                  //
////////////////////////////////////////////////////////////////////////////

void readInFormedParticles(int Npart, float *r0, float *r1, float *r2, float *r3, float *p0, float *p1, float *p2, float *p3, float *mi, float *gi, float *bi){
    
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
    
    FILE *allsetfile;
    char allsetname[255];
    sprintf(allsetname, "%s.dat", "output/AllSets");
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
    
}

