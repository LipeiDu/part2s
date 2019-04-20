#include <H5Cpp.h>
#include <H5File.h>

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

void writeSourcesASCII(int n, int Nx, int Ny, int Nn, int Ntot, float tau, float dt, float dx, float dy, float dn, float *Sb, float *St, float *Sx, float *Sy, float *Sn, float *Norm){
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
                
                float x   = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
                float y   = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
                float eta = ((float)k - ((float)Nn - 1.0)/2.0) * dn;
                
                int s = i + j * (Nx) + k * (Nx * Ny);
                
                fprintf(sourcefile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, St[s], Sx[s], Sy[s], Sn[s], Sb[s]);
                
                *Norm = *Norm + Sb[s] * tau;
            } // for (int k )
        } //for (int j)
    } //for (int i )
    
    printf("Numerical Nb = %lf\n", (*Norm) * dt * dx * dy * dn);
    fclose(sourcefile);
}

void writeTensorsASCII(int n, int Nx, int Ny, int Nn, int Ntot, float tau, float dt, float dx, float dy, float dn, float *Sb, float *St, float *Sx, float *Sy, float *Sn, float *Norm){
#ifdef INITIAL_TENSOR
    FILE *tmunufile;
    char fjname[255];
    sprintf(fjname, "%s%d.dat", "output/Tmunu", n);
    tmunufile = fopen(fjname, "w");
    
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
                
                fprintf(tmunufile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", x, y, eta, energyDensity[s], flowVelocity[0][s], flowVelocity[1][s], flowVelocity[2][s], flowVelocity[3][s]);
                
            }
        }
    }
    
    fclose(tmunufile);
#endif
}

