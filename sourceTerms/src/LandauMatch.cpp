
////////////////////////////////////////////////////////////////////////////
//                            Landau Matching                             //
////////////////////////////////////////////////////////////////////////////

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#define REGULATE 1 // 1 to regulate flow in dilute regions
#define GAMMAMAX 100.0
#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

double equilibriumPressure(double e);

void solveEigenSystem(float **stressTensor, float *energyDensity, float **flowVelocity, float *pressure, float TAU, int DIM_X, int DIM_Y, int DIM_ETA, int DIM, float DX, float DY)
{
    //********************************************************
    // Solve eigenvalue problem in the Landau Frame
    //********************************************************
    
    for (int is = 0; is < DIM; is++)
    {
        gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
        
        //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
        Tmunu = gsl_matrix_alloc(4,4);
        gsl_matrix *gmunu;
        gmunu = gsl_matrix_alloc(4,4);
        gsl_matrix_complex *eigen_vectors;
        eigen_vectors = gsl_matrix_complex_alloc(4,4);
        gsl_vector_complex *eigen_values;
        eigen_values = gsl_vector_complex_alloc(4);
        
        //set the values of the energy momentum tensor
        gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is]); //tau,tau
        gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][is]); //tau,x
        gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][is]); //tau,y
        gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][is]); //tau,eta
        gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][is]); //x,x
        gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][is]); //x,y
        gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][is]); //x,eta
        gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][is]); //y,y
        gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][is]); //y,eta
        gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][is]); //eta,eta
        gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][is]); //x,tau
        gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][is]); //y,tau
        gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][is]); //eta,tau
        gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][is]); //y,x
        gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][is]); //eta,x
        gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][is]); //eta,y
        
        //set the values of the "metric"; not really the metric, but the numerical constants
        //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
        //note factors of TAU appropriate for milne coordinates g_(mu.nu) = diag(1,-1,-1,-TAU^2)
        gsl_matrix_set(gmunu, 0, 0, 1.0); //tau,tau
        gsl_matrix_set(gmunu, 0, 1, -1.0); //tau,x
        gsl_matrix_set(gmunu, 0, 2, -1.0); //tau,y
        gsl_matrix_set(gmunu, 0, 3, -1.0*TAU*TAU); //tau,eta
        gsl_matrix_set(gmunu, 1, 0, 1.0); //x,tau
        gsl_matrix_set(gmunu, 1, 1, -1.0); //x,x
        gsl_matrix_set(gmunu, 1, 2, -1.0); //x,y
        gsl_matrix_set(gmunu, 1, 3, -1.0*TAU*TAU); //x,eta
        gsl_matrix_set(gmunu, 2, 0, 1.0); //y,tau
        gsl_matrix_set(gmunu, 2, 1, -1.0); //y,x
        gsl_matrix_set(gmunu, 2, 2, -1.0); //y,y
        gsl_matrix_set(gmunu, 2, 3, -1.0*TAU*TAU); //y,eta
        gsl_matrix_set(gmunu, 3, 0, 1.0); //eta,tau
        gsl_matrix_set(gmunu, 3, 1, -1.0); //eta,x
        gsl_matrix_set(gmunu, 3, 2, -1.0); //eta,y
        gsl_matrix_set(gmunu, 3, 3, -1.0*TAU*TAU); //eta,eta
        
        //lower one index of the stress tensor; save it to the same matrix to save memory
        gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
        gsl_eigen_nonsymmv_workspace *eigen_workspace;
        eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
        gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
        gsl_eigen_nonsymmv_free(eigen_workspace);
        
        //***does this have a solution for energy density and flow at every point?
        int eigenvalue_exists = 0;
        for (int i = 0; i < 4; i++)
        {
            gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);
            
            if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
            {
                gsl_complex v0 = gsl_matrix_complex_get(eigen_vectors, 0 , i);
                gsl_complex v1 = gsl_matrix_complex_get(eigen_vectors, 1 , i);
                gsl_complex v2 = gsl_matrix_complex_get(eigen_vectors, 2 , i);
                gsl_complex v3 = gsl_matrix_complex_get(eigen_vectors, 3 , i);
                
                if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
                {
                    double minkowskiLength = GSL_REAL(v0)*GSL_REAL(v0) - (GSL_REAL(v1)*GSL_REAL(v1) + GSL_REAL(v2)*GSL_REAL(v2) + TAU*TAU*GSL_REAL(v3)*GSL_REAL(v3));
                    double factor = 1.0 / sqrt(minkowskiLength);
                    
                    if (GSL_REAL(v0) < 0) factor=-factor;
                    
                    //ignore eigenvectors with gamma >~ 60
                    if ( (GSL_REAL(v0) * factor) < GAMMAMAX)
                    {
                        eigenvalue_exists = 1;
                        energyDensity[is] = GSL_REAL(eigenvalue);
                        flowVelocity[0][is] = GSL_REAL(v0) * factor;
                        flowVelocity[1][is] = GSL_REAL(v1) * factor;
                        flowVelocity[2][is] = GSL_REAL(v2) * factor;
                        flowVelocity[3][is] = GSL_REAL(v3) * factor;
                    }
                    
                } // if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
            } // if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
        } //for (int i = 0; i < 4; ...)
        
        if (eigenvalue_exists == 0)
        {
            //in dilute regions where we can't find a timelike eigenvector, set e = 0, u^t = 1, u^x=u^y=u^n=0
            energyDensity[is] = 0.0;
            flowVelocity[0][is] = 1.0;
            flowVelocity[1][is] = 0.0;
            flowVelocity[2][is] = 0.0;
            flowVelocity[3][is] = 0.0;
        }
    } // for (int is; is < DIM; ...)
    
    
    //********************************************************
    // regulation
    //********************************************************
    
    //now regulate the flow velocity by a smoothing procedure. Flow can be too large in dilute regions, cause hydro to crash...
    //this method doesnt yield a smooth profile
    
    //try scaling the flow velocity by a smooth profile which goes to zero after some finite radius
    if (REGULATE)
    {
        printf("Regulating flow velocity profile in dilute regions \n");
        for (int is = 0; is < DIM; is++){
            
            int i = is / (DIM_Y * DIM_ETA);
            int j = (is - (DIM_Y * DIM_ETA * i))/ DIM_ETA;
                    
            double x = (i - (DIM_X-1)/2.)*DX;
            double y = (j - (DIM_Y-1)/2.)*DY;
        
            double r = sqrt(x*x + y*y);
        
            float R_WIDTH = 0.6;
            float R_FLAT = 4.5;
            float arg = (-1.0) * (r - R_FLAT) * (r - R_FLAT) / (2.0 * R_WIDTH * R_WIDTH);
            arg = arg * THETA_FUNCTION(r - R_FLAT);
        
            flowVelocity[1][is] = flowVelocity[1][is] * exp(arg);
            flowVelocity[2][is] = flowVelocity[2][is] * exp(arg);
        
            flowVelocity[0][is] = sqrt( 1 + flowVelocity[1][is]*flowVelocity[1][is] + flowVelocity[2][is]*flowVelocity[2][is] + TAU*TAU*flowVelocity[3][is]*flowVelocity[3][is]);
        }
    }
    
    //********************************************************
    // Initialize pressure
    //********************************************************
    
    for (int is = 0; is < DIM; is++){
                pressure[is] = equilibriumPressure(energyDensity[is]);
    }
}

// PI = -1/3 * (T^(mu)_(mu) - epsilon) - p
// T^(mu)_(mu) = T^(0,0) - T^(1,1) - T^(2,2) - (TAU^2)T^(3,3)
void calculateBulkPressure(float **stressTensor, float *energyDensity, float *pressure, float *bulkPressure, int DIM, float TAU)
{
    for (int is = 0; is < DIM; is++)
    {
        float a =  stressTensor[0][is] - stressTensor[4][is] - stressTensor[7][is] - TAU*TAU*stressTensor[9][is];
        bulkPressure[is] = (-1.0/3.0) * (a - energyDensity[is]) - pressure[is];
    }
}

// pi^(mu,nu) = T^(mu,nu) - epsilon * u^(mu)u^(nu) + (P + PI) * (g^(mu,nu) - u^(mu)u^(nu))
void calculateShearViscTensor(float **stressTensor, float *energyDensity, float **flowVelocity, float *pressure, float *bulkPressure, float **shearTensor, int DIM, float TAU)
{
    for (int is = 0; is < DIM; is++)
    {
        //calculate ten components : upper triangular part
        float b = energyDensity[is] + pressure[is] + bulkPressure[is];
        float c = pressure[is] + bulkPressure[is];
        shearTensor[0][is] = stressTensor[0][is] - flowVelocity[0][is] * flowVelocity[0][is] * b + c; //pi^(tau,tau)
        shearTensor[1][is] = stressTensor[1][is] - flowVelocity[0][is] * flowVelocity[1][is] * b; //pi^(tau,x)
        shearTensor[2][is] = stressTensor[2][is] - flowVelocity[0][is] * flowVelocity[2][is] * b; //pi^(tau,y)
        shearTensor[3][is] = stressTensor[3][is] - flowVelocity[0][is] * flowVelocity[3][is] * b; //pi^(tau,eta)
        shearTensor[4][is] = stressTensor[4][is] - flowVelocity[1][is] * flowVelocity[1][is] * b - c; //pi^(x,x)
        shearTensor[5][is] = stressTensor[5][is] - flowVelocity[1][is] * flowVelocity[2][is] * b; //pi^(x,y)
        shearTensor[6][is] = stressTensor[6][is] - flowVelocity[1][is] * flowVelocity[3][is] * b; //pi^(x,eta)
        shearTensor[7][is] = stressTensor[7][is] - flowVelocity[2][is] * flowVelocity[2][is] * b - c; //pi^(y,y)
        shearTensor[8][is] = stressTensor[8][is] - flowVelocity[2][is] * flowVelocity[3][is] * b; //pi^(y,eta)
        shearTensor[9][is] = stressTensor[9][is] - flowVelocity[3][is] * flowVelocity[3][is] * b - c * (1.0/(TAU*TAU)); //pi^(eta,eta)
    }
}

// n_B = u^(mu)j_(mu)
void calculateBaryonDensity(float *baryonDensity, float **baryonCurrent, float **flowVelocity, int DIM, float TAU)
{
    for (int is = 0; is < DIM; is++)
    {
        baryonDensity[is] = (flowVelocity[0][is] * baryonCurrent[0][is]) - (flowVelocity[1][is] * baryonCurrent[1][is]) - (flowVelocity[2][is] * baryonCurrent[2][is]) - (TAU * TAU * flowVelocity[3][is] * baryonCurrent[3][is]);
    }
}

// V^(mu) = j^(mu) - n_B * u^(mu)
void calculateBaryonDiffusion(float **baryonDiffusion, float **baryonCurrent, float *baryonDensity, float **flowVelocity, int DIM)
{
    for (int ivar = 0; ivar < 4; ivar++)
    {
        for (int is = 0; is < DIM; is++)
        {
            baryonDiffusion[ivar][is] = baryonCurrent[ivar][is] - (baryonDensity[is] * flowVelocity[ivar][is]);
        }
    }
}

// Equation of State
double equilibriumPressure(double e) {
    // Equation of state from the Wuppertal-Budapest collaboration
    double e1 = (double)e;
    double e2 = e*e;
    double e3 = e2*e;
    double e4 = e3*e;
    double e5 = e4*e;
    double e6 = e5*e;
    double e7 = e6*e;
    double e8 = e7*e;
    double e9 = e8*e;
    double e10 = e9*e;
    double e11 = e10*e;
    double e12 = e11*e;
    
    double a0 = -0.25181736420168666;
    double a1 = 9737.845799644809;
    double a2 = 1.077580993288114e6;
    double a3 = 3.1729694865420084e6;
    double a4 = 1.6357487344679043e6;
    double a5 = 334334.4309240126;
    double a6 = 41913.439282708554;
    double a7 = 6340.448389300905;
    double a8 = 141.5073484468774;
    double a9 = 0.7158279081255019;
    double a10 = 0.0009417586777847889;
    double a11 = 3.1188455176941583e-7;
    double a12 = 1.9531729608963267e-11;
    double a = (double)fma(a12,e12,fma(a11,e11,fma(a10,e10,fma(a9,e9,fma(a8,e8,fma(a7,e7,fma(a6,e6,fma(a5,e5,fma(a4,e4,fma(a3,e3,fma(a2,e2,fma(a1,e1,a0))))))))))));
    
    double b0 = 45829.44617893836;
    double b1 = 4.0574329080826794e6;
    double b2 = 2.0931169138134286e7;
    double b3 = 1.3512402226067686e7;
    double b4 = 1.7851642641834426e6;
    double b5 = 278581.2989342773;
    double b6 = 26452.34905933697;
    double b7 = 499.04919730607065;
    double b8 = 2.3405487982094204;
    double b9 = 0.002962497695527404;
    double b10 = 9.601103399348206e-7;
    double b11 = 5.928138360995685e-11;
    double b12 = 3.2581066229887368e-18;
    double b = (double)fma(b12,e12,fma(b11,e11,fma(b10,e10,fma(b9,e9,fma(b8,e8,fma(b7,e7,fma(b6,e6,fma(b5,e5,fma(b4,e4,fma(b3,e3,fma(b2,e2,fma(b1,e1,b0))))))))))));
    return a/b;
}


