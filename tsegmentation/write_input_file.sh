#!/bin/bash

filename=$1
Kn=$2 ; #default 5e4
Kt=$3 ; #default 5e4
Bn=$4 ; #default 2e5
Bt=$5 ; #default 2e5
Bm=$6 ; #default 2e5
Eps=$7 ; #default 0.01
maxt=$8 ; #default 1.0
cat << EOF > $filename
0.05    = Verlet:      Verlet distance for optimization
voronoi = ptype:       Particle type (voronoi, sphere,rice and tetrahedra mesh so far)
flexion = test:        Type of test (tensile or vibration)
2       = RenderVideo: 1 if video should be render, 0 otherwise
${Kn}   = Kn:          Normal stiffness
${Kt}   = Kt:          Tangential stiffness
-0.2    = Gn:          Normal dissipative coefficient (contacts)
-0.0    = Gt:         Tangential dissipative coefficient (contacts)
0.3     = Mu:          Microscopic friction coefficient
${Bn}   = Bn:          Cohesion Normal stiffness
${Bn}   = Bt:          Cohesion Tangential stiffness
${Bm}   = Bm:          Cohesion Torque stiffness
${Eps}  = Eps:         Cohesion Threshold (if negative, the bonds will never break)
0.1     = R:           Spheroradius
1000    = seed:        Seed of the random generator
5.0e-5  = dt:          Time step
0.1     = dtOut:       Time step for output
30.0    = Tf:          Final time for the test
33.0    = Lx:          Lx (length of sphere cube)
11.0    = Ly:          Ly (irrelevant for spheres)
23.0    = Lz:          Lz (irrelevant for spheres)
1.0     = dx:          dx (Spacing for brick geometry)
10      = nx:          nx (number of spheres per cube size)
3       = ny:          ny (irrelevant for spheres)
7       = nz:          nz (irrelevant for spheres)
1.5     = Rc:          Radius for the cylinders
3.0     = rho:         rho
-600.0  = Am:          vibration force amplitude
0.3     = ome:         Frequency of vibration
-0.105   = ex:          Velocity of the compressing plane (positive for extension, negative for compression)
1.0e-3  = SEFthr:      Strain energy field threshold for breaking particles
0       = Restart:     Restart state, 0 Don't restart, 1, 2, ... restart number
${maxt} = max_time:    Maximum time at which compression should run, after this the system equilibrates
200     = idx_init:    Last index from previous simulation
1       = cohesion:    Whether or not to add cohesion between particles
EOF
