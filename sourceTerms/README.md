# part2s

## code structer

### For each event from UrQMD:
* Get the particle list in Cartisian coordinates
* Calculate the Milne coordinates at formation time, write in Milne.dat file
* Calculate the Center of Mass of the event at the same time
* Get the particle list in Cartisian coordinates again, and calculate the participant plane angle with respect to the CM
* Open the particle list in Milne.dat file, recenter and rotate x, y and px, py with CM and participant plane angle
* Write the new coordinates in a file (all following events in this file to be averaged over for source terms)
### With the particle list from all events
* Read in the coordinates and store in array
* Calculate source terms


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
    float kernel = 1.0/(ch * ch) * numerator; // [1]
    
    float tauInv = 1.0 / tau;
    
    // pi[m][4] is [GeV] by defination above
    Sb = Sb + kernel * b_i; // [1]
    St = St + kernel * p0_d[m]; // [GeV]
    Sx = Sx + kernel * p1_d[m]; // [GeV]
    Sy = Sy + kernel * p2_d[m]; // [GeV]
    Sn = Sn + kernel * p3_d[m] * tauInv; // [GeV/fm] caution, definition and tau here
    
    float facN = prefactor * nevInv;
    float facHN = prefactor * hbarcNevInv; //doing average and unit conversion
    
    // dtau*dx*dy*tau*deta
    Sb_d[tid] = facN  * Sb * tauInv; // [1/fm^4] = [1/fm^3] * [1] * [1/fm]
    St_d[tid] = facHN * St * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
    Sx_d[tid] = facHN * Sx * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
    Sy_d[tid] = facHN * Sy * tauInv; // [1/fm^5] = [1/(fm^4*GeV)] * [GeV] * [1/fm]
    Sn_d[tid] = facHN * Sn * tauInv; // [1/fm^6] = [1/(fm^4*GeV)] * [GeV/fm] * [1/fm]
```
