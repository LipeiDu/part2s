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

void writeSourcesHDF5(int n, int Nx, int Ny, int Nn, int Ntot, float *Sall, float *Sb, float *St, float *Sx, float *Sy, float *Sn){
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
    char source_fname[255];
    sprintf(source_fname, "%s%d.h5", "output/Sources", n);
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
}

void writeSourcesASCII(int n, int Nx, int Ny, int Nn, int Ntot, float tau, float dt, float dx, float dy, float dn, float *Sb, float *St, float *Sx, float *Sy, float *Sn, float *Sbtotal, float *S0total, int nev){
    FILE *sourcefile;
    char finame[255];
    
    sprintf(finame, "%s%d.dat", "output/Sources", n);
    sourcefile = fopen(finame, "w");
    
    float dV = dt * dx * dy * dn;
    
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
                
                *S0total += tau * StTot(tau, eta, St[s], Sx[s], Sy[s], Sn[s]) * dV;
                *Sbtotal += tau * Sb[s] * dV;
                
            } // for (int k )
        } //for (int j)
    } //for (int i )
    
    *S0total *= nev * hbarc;
    *Sbtotal *= nev;
    
    printf("Sbtotal = %lf,\t S0total = %lf\n", *Sbtotal, *S0total);
    
    fclose(sourcefile);
}

void writeTensorsASCII(int n, int Nx, int Ny, int Nn, int Ntot, float tau, float dt, float dx, float dy, float dn){
#ifdef INITIAL_TENSOR
    FILE *tmunufile, *pifile, *baryonfile;
    char finame[255], fjname[255], fkname[255];
    
    sprintf(finame, "%s%d.dat", "output/BaryonVariables", n);
    baryonfile = fopen(finame, "w");
    sprintf(fjname, "%s%d.dat", "output/PrimaryVariables", n);
    tmunufile = fopen(fjname, "w");
    sprintf(fkname, "%s%d.dat", "output/DissipativeTerms", n);
    pifile = fopen(fkname, "w");
    
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
                
                // x,y,eta,n,Nt,Nx,nt,nx
                fprintf(baryonfile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, baryonDensity[s], baryonCurrent[0][s], baryonCurrent[1][s], baryonDiffusion[0][s], baryonDiffusion[1][s]);
                
                // x,y,eta,e,ut,ux,uy,un
                fprintf(tmunufile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, energyDensity[s], flowVelocity[0][s], flowVelocity[1][s], flowVelocity[2][s], flowVelocity[3][s]);
                
                // x,y,eta,pitt,pitx,pixy,pinn,Pi
                fprintf(pifile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, shearTensor[0][s], shearTensor[1][s], shearTensor[5][s], shearTensor[9][s], bulkPressure[s]);

            }
        }
    }
    
    fclose(baryonfile);
    fclose(tmunufile);
    fclose(pifile);
#endif
}

