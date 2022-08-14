#!/bin/bash

for kn in 1e-1 1e0 1e1 1e2 1e3 1e4 1e7
do
	echo "Running kn=${kn}"
	time ./simple_breaking_restart_equilv1 brick_restart_kn${kn}_equil 10 > simple_breaking_restart_kn${kn}_equil.log
	./generate_cylinder_logs.sh simple_breaking_restart_kn${kn}_equil.log 1.0 1.0_kn${kn}
	./move_files.sh flexion_1.0s_gm1.0_kn${kn}_equal_radii
done
