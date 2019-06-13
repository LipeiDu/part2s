
#include "ParameterReader.cpp"

// refered to Duke-QCD Trento code

void calculateEccentricity(float *energyDensity, int n, int Nx, int Ny, int Nn, float dx, float dy, float gmax, float sigma){
    
    //****************************************************
    // find the center of mass
    //****************************************************
    
    double eTot = 0.;
    double eTotX = 0.;
    double eTotY = 0.;
    
    // here we calculate eccentricities time step by time step, the factor tau*dt*dn or tau*dt*dx*dy*dn is a common factor, will cancel out
    
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nn; ++k)
            {
                
                float x   = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
                float y   = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
                
                int s = i + j * (Nx) + k * (Nx * Ny);

                eTot += energyDensity[s];

                eTotX += energyDensity[s] * x;
                eTotY += energyDensity[s] * y;
            }
        }
    }
    
    double xCM = 0.;
    double yCM = 0.;
    
    xCM = eTotX / eTot;
    yCM = eTotY / eTot;
    
    
    //****************************************************
    // calculate the eccentricities, see DUKE-QCD TRENTO
    //****************************************************
    
    // real part, imaginary part and weight
    float e2r = 0.0;
    float e2i = 0.0;
    float e2w = 0.0;
    float e3r = 0.0;
    float e3i = 0.0;
    float e3w = 0.0;
    float e4r = 0.0;
    float e4i = 0.0;
    float e4w = 0.0;
    float e5r = 0.0;
    float e5i = 0.0;
    float e5w = 0.0;

    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nn; ++k)
            {
                
                float xr   = ((float)i - ((float)Nx - 1.0)/2.0) * dx;
                float yr   = ((float)j - ((float)Ny - 1.0)/2.0) * dy;
                
                int s = i + j * (Nx) + k * (Nx * Ny);
                
                float e = energyDensity[s];
                if (e < 1.e-10)
                    continue;
                
                float x  = xr - xCM;
                float x2 = x  * x;
                float x3 = x2 * x;
                float x4 = x2 * x2;
                
                float y  = yr - yCM;
                float y2 = y  * y;
                float y3 = y2 * y;
                float y4 = y2 * y2;
                
                float r2 = x2 + y2;
                float r  = sqrt(r2);
                float r4 = r2 * r2;
                
                float xy = x * y;
                float x2y2 = x2 * y2;
                
                e2r += e * (y2 - x2);
                e2i += e * 2.0 *xy;
                e2w += e * r2;
                
                e3r += e * (y3 - 3.0 * y * x2);
                e3i += e * (3.0 * x * y2 - x3);
                e3w += e * r2 * r;
                
                e4r += e * (x4 + y4 - 6.0 * x2y2);
                e4i += e * 4.0 * xy * (y2 - x2);
                e4w += e * r4;
                
                e5r += e * y * (5.0 * x4 - 10.0 * x2y2 + y4);
                e5i += e * x * (x4 - 10.0 * x2y2 + 5.0 * y4);
                e5w += e * r4 * r;
            }
        }
    }
    
    float eccentricity2 = sqrt(e2r * e2r + e2i * e2i) / e2w;
    float eccentricity3 = sqrt(e3r * e3r + e3i * e3i) / e3w;
    float eccentricity4 = sqrt(e4r * e4r + e4i * e4i) / e4w;
    float eccentricity5 = sqrt(e5r * e5r + e5i * e5i) / e5w;
    
    //****************************************************
    // write the file
    //****************************************************

    FILE *eccfile;
    char filename[255];
    
    sprintf(filename, "%s%d.dat", "output/eccentricities", n);
    eccfile = fopen(filename, "w");

    fprintf(eccfile, "%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", gmax, sigma, eccentricity2, eccentricity3, eccentricity4, eccentricity5);
    
    fclose(eccfile);
}

