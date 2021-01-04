#include <H5Cpp.h>
#include <H5File.h>
#include "ParameterReader.cpp"


float TttTot(float tau, float eta, float Ttt, float Ttx, float Tty, float Ttn, float Txx, float Txy, float Txn, float Tyy, float Tyn, float Tnn){
    
    float cosheta = cosh(eta);
    float sinheta = sinh(eta);
    
    return cosheta * cosheta * Ttt + 2 * tau * cosheta * sinheta * Ttn + tau * tau * sinheta * sinheta * Tnn;
}

float StTot(float tau, float eta, float St, float Sx, float Sy, float Sn){
    
    float cosheta = cosh(eta);
    float sinheta = sinh(eta);
    
    return cosheta * St + tau * sinheta * Sn;
}

////////////////////////////////////////////////////////////////////////////
// HDF5 files
////////////////////////////////////////////////////////////////////////////

void writeSourcesHDF5(int n, int Nx, int Ny, int Nn, int Ntot, float *Sall, float *Sb, float *St, float *Sx, float *Sy, float *Sn){
    
    int N2 = 5;
    
    //compress all data into the 1d array to pass to hdf5 writer
    for (int is = 0; is < Ntot; is++)
    {
   
        float Sarray[5] = {St[is], Sx[is], Sy[is], Sn[is], Sb[is]};
        
        for (int i = 0; i < N2; i++){
            
            Sall[N2 * is + i] = Sarray[i];
            
        }
    }
    
    //printf("Writing source terms to file...\n\n");
    char source_fname[255];
    sprintf(source_fname, "%s%d.h5", "output/Sources", n);
    H5::H5File file(source_fname, H5F_ACC_TRUNC);
    
    hsize_t dimsf[2];
    dimsf[0] = Nx * Ny * Nn;
    dimsf[1] = N2;
    
    H5::DataSpace dataspace(2, dimsf);
    
    H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    dataset.write(Sall, H5::PredType::NATIVE_FLOAT);
    
}

//Tmunu
void writeStressHDF5(int n, int Nx, int Ny, int Nn, int Ntot, float *stressAll, float **stressTensor){
    
    int N2 = 10;
    
    //compress all data into the 1d array to pass to hdf5 writer
    for (int is = 0; is < Ntot; is++)
    {
        float Tarray[10] = {stressTensor[0][is], stressTensor[1][is], stressTensor[2][is], stressTensor[3][is], stressTensor[4][is], stressTensor[5][is], stressTensor[6][is], stressTensor[7][is], stressTensor[8][is], stressTensor[9][is]};
        
        for (int i = 0; i < N2; i++){
            
            stressAll[N2 * is + i] = Tarray[i];
            
        }
    }
    
    //printf("Writing source terms to file...\n\n");
    char source_fname[255];
    sprintf(source_fname, "%s%d.h5", "output/stressTensor", n);
    H5::H5File file(source_fname, H5F_ACC_TRUNC);
    
    hsize_t dimsf[2];
    dimsf[0] = Nx * Ny * Nn;
    dimsf[1] = N2;
    
    H5::DataSpace dataspace(2, dimsf);
    
    H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    dataset.write(stressAll, H5::PredType::NATIVE_FLOAT);
    
}

//shear stress
void writeShearHDF5(int n, int Nx, int Ny, int Nn, int Ntot, float *shearAll, float **shearTensor){
        
    int N2 = 10;
    
    //compress all data into the 1d array to pass to hdf5 writer
    for (int is = 0; is < Ntot; is++)
    {
        
        float Sharray[10] = {shearTensor[0][is], shearTensor[1][is], shearTensor[2][is], shearTensor[3][is], shearTensor[4][is], shearTensor[5][is], shearTensor[6][is], shearTensor[7][is], shearTensor[8][is], shearTensor[9][is]};
        
        for (int i = 0; i < N2; i++){
            
            shearAll[N2 * is + i] = Sharray[i];
            
        }
    }
    
    //printf("Writing source terms to file...\n\n");
    char source_fname[255];
    sprintf(source_fname, "%s%d.h5", "output/shearTensor", n);
    H5::H5File file(source_fname, H5F_ACC_TRUNC);
    
    hsize_t dimsf[2];
    dimsf[0] = Nx * Ny * Nn;
    dimsf[1] = N2;
    
    H5::DataSpace dataspace(2, dimsf);
    
    H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    dataset.write(shearAll, H5::PredType::NATIVE_FLOAT);   
    
}

//energy, pressure, temperature, bulk, flow velocity
void writePrimaryHDF5(int n, int Nx, int Ny, int Nn, int Ntot, float *primaryAll, float *energyDensity, float *pressure, float *temperature, float *bulkPressure, float **flowVelocity){
    
    int N2 = 8;
    
    //compress all data into the 1d array to pass to hdf5 writer
    for (int is = 0; is < Ntot; is++)
    {
        
        float Parray[8] = {energyDensity[is], pressure[is], temperature[is], bulkPressure[is], flowVelocity[0][is], flowVelocity[1][is], flowVelocity[2][is], flowVelocity[3][is]};
        
        for (int i = 0; i < N2; i++){
            
            primaryAll[N2 * is + i] = Parray[i];
            
        }
    }
    
    //printf("Writing source terms to file...\n\n");
    char source_fname[255];
    sprintf(source_fname, "%s%d.h5", "output/primaryVariables", n);
    H5::H5File file(source_fname, H5F_ACC_TRUNC);
    
    hsize_t dimsf[2];
    dimsf[0] = Nx * Ny * Nn;
    dimsf[1] = N2;
    
    H5::DataSpace dataspace(2, dimsf);
    
    H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    dataset.write(primaryAll, H5::PredType::NATIVE_FLOAT);  
    
}

//baryon density, net current, diffusion
void writeBaryonHDF5(int n, int Nx, int Ny, int Nn, int Ntot, float *baryonAll, float *baryonDensity, float **baryonCurrent, float **baryonDiffusion){
    
    int N2 = 9;
    
    //compress all data into the 1d array to pass to hdf5 writer
    for (int is = 0; is < Ntot; is++)
    {
        
        float Barray[9] = {baryonDensity[is], baryonCurrent[0][is], baryonCurrent[1][is], baryonCurrent[2][is], baryonCurrent[3][is], baryonDiffusion[0][is], baryonDiffusion[1][is], baryonDiffusion[2][is], baryonDiffusion[3][is]};
        
        for (int i = 0; i < N2; i++){
            
            baryonAll[N2 * is + i] = Barray[i];
            
        }
    }
    
    //printf("Writing source terms to file...\n\n");
    char source_fname[255];
    sprintf(source_fname, "%s%d.h5", "output/baryonVariables", n);
    H5::H5File file(source_fname, H5F_ACC_TRUNC);
    
    hsize_t dimsf[2];
    dimsf[0] = Nx * Ny * Nn;
    dimsf[1] = N2;
    
    H5::DataSpace dataspace(2, dimsf);
    
    H5::DataType datatype(H5::PredType::NATIVE_FLOAT);
    H5::DataSet dataset = file.createDataSet("data", datatype, dataspace);
    dataset.write(baryonAll, H5::PredType::NATIVE_FLOAT);
    
}

////////////////////////////////////////////////////////////////////////////
// dat files
////////////////////////////////////////////////////////////////////////////

void writeSourcesASCII(int n, int Nx, int Ny, int Nn, int Ntot, float tau, float dt, float dx, float dy, float dn, float *Sb, float *St, float *Sx, float *Sy, float *Sn, float *Sbn, float *Stn, int nev){
    //FILE *sourcefile;
    //char finame[255];
    
//     sprintf(finame, "%s%d.dat", "output/Sources", n);
//     sourcefile = fopen(finame, "w");
    
    float dV = dt * dx * dy * dn;
    
    for (int k = 0; k < Nn; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                
                float x   = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
                float y   = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
                float eta = ((float)k - ((float)Nn - 1.0)/2.0) * dn;
                
                int s = i + j * (Nx) + k * (Nx * Ny);
                
                //fprintf(sourcefile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, St[s], Sx[s], Sy[s], Sn[s], Sb[s]);
                
                //if(isinf(Sb[s]))
                //    printf("Sb[s] is inf,\t xyn = %lf\t%lf\t%lf\n", Sb[s], x,y,eta);
                
                //printf("Before: Sbtotal is = %lf\n", *Sbtotal);
                
                *Stn += tau * StTot(tau, eta, St[s], Sx[s], Sy[s], Sn[s]) * dV;
                *Sbn += tau * Sb[s] * dV;
                
                
                
                //if(isinf(*Sbtotal))
                //    printf("Sbtotal is inf,\t Sb[s],tau,dv = %lf\t%lf\t%lf\n", Sb[s], tau, dV);
                
            } // for (int k )
        } //for (int j)
    } //for (int i )
    
    printf("Sbn is = %lf, Stn is = %lf\n", *Sbn, *Stn);
    
    //*S0total *= nev * hbarc;
    //*Sbtotal *= nev;
    
    //printf("Sbtotal = %lf,\t S0total = %lf\n", *Sbtotal, *S0total);
    
    //fclose(sourcefile);
}

void writeTensorsASCII(int n, int Nx, int Ny, int Nn, int Ntot, float tau, float dt, float dx, float dy, float dn){
#ifdef INITIAL_TENSOR
    FILE *tmunufile, *pifile, *baryonfile, *tensorfile;
    char finame[255], fjname[255], fkname[255], flname[255];
    
    sprintf(finame, "%s%d.dat", "output/BaryonVariables", n);
    baryonfile = fopen(finame, "w");
    sprintf(fjname, "%s%d.dat", "output/PrimaryVariables", n);
    tmunufile = fopen(fjname, "w");
    sprintf(fkname, "%s%d.dat", "output/DissipativeTerms", n);
    pifile = fopen(fkname, "w");
    sprintf(flname, "%s%d.dat", "output/Tmunu", n);
    tensorfile = fopen(flname, "w");
    
    for (int k = 0; k < Nn; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                
                float x   = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
                float y   = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
                float eta = ((float)k - ((float)Nn - 1.0)/2.0) * dn;
                
                int s = i + j * (Nx) + k * (Nx * Ny);
                
                // x,y,eta,n,Nt,Nx,nt,nx
                fprintf(baryonfile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, baryonDensity[s], baryonCurrent[0][s], baryonCurrent[1][s], baryonDiffusion[0][s], baryonDiffusion[1][s]);
                
                // x,y,eta,e,ut,ux,uy,un
                fprintf(tmunufile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, energyDensity[s], flowVelocity[0][s], flowVelocity[1][s], flowVelocity[2][s], flowVelocity[3][s]);
                
                // x,y,eta,pitt,pitx,pixy,pinn,Pi
                fprintf(pifile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, shearTensor[0][s], shearTensor[1][s], shearTensor[5][s], shearTensor[9][s], bulkPressure[s]);
                
                // x,y,eta,...
                fprintf(tensorfile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, energyDensity[s], pressure[s], flowVelocity[0][s], flowVelocity[1][s], flowVelocity[2][s], flowVelocity[3][s], shearTensor[0][s], shearTensor[1][s], shearTensor[2][s], shearTensor[3][s], shearTensor[4][s], shearTensor[5][s], shearTensor[6][s], shearTensor[7][s], shearTensor[8][s], shearTensor[9][s], bulkPressure[s]);

            }
        }
    }
    
    fclose(baryonfile);
    fclose(tmunufile);
    fclose(pifile);
    fclose(tensorfile);
#endif
}

