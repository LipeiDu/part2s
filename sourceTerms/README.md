# part2s

## notation and demension

Smearing kernel
```
    float SigInv = 1.0/(2.0*sigma*sigma);
    float SignInv = 1.0/(2.0*sigman*sigman);
    float d_tauInv = 1.0/delta_tau;
    float nevInv = 1.0/(float)nev;
    float hbarcNevInv = nevInv/hbarc;

    float prefac = 1.0/(2.0 * (2.0*PI*sigma*sigma) * sqrt(2.0*PI*sigman*sigman));
    float prefactor = d_tauInv * prefac; 
    
    float distt = fabs(rr[0]-r0_d[m]);
    float ddx = fabs(rr[1]-r1_d[m]);
    float ddy = fabs(rr[2]-r2_d[m]);
    float disttrs = ddx*ddx + ddy*ddy;
    float disttr = sqrt(disttrs);
    float distn = fabs(rr[3]-r3_d[m]);
    float dist = -(disttrs * SigInv + distn*distn * SignInv);
    
    float numerator = exp(dist);
    float delta = distt * d_tauInv;
    float ch = cosh(delta);
    float kernel = 1.0/(ch * ch) * numerator;
    
    float tauInv = 1.0 / tau;
    
    Sb = Sb + kernel * b_i;
    St = St + kernel * p0_d[m];
    
    float facN = prefactor * nevInv;
    float facHN = prefactor * hbarcNevInv; //doing average and unit conversion
    
    Sb_d[tid] = facN  * Sb * tauInv;
    St_d[tid] = facHN * St * tauInv;
```
