#!/bin/bash

#run_0000
#arguments are inp_file Ncores RestartIndex paticletype
#paramsline="150s_maxt3.0s_gm4.0_kn5e2_bn1e4_1e-4_v0.025"
#inpfile="brick_restart_iter"
#breaktime="150.000000"
paramsline="30s_maxt1.0s_gm4.0_kn5e4_bn2e5_1e-4_v0.105"
inpfile="brick_restart_fast"
breaktime="29.900100"
#"1.000100"
#"29.900100"
time ./simple_breaking_restart_equilv2 ${inpfile} 10 0 "voronoi" > simple_breaking_restart_${paramsline}_iter0.log
./generate_cylinder_logs.sh simple_breaking_restart_${paramsline}_iter0.log ${paramsline}
./move_files.sh run_0000_${paramsline}
./cut_particles brick_restart_geometry 10 0 run_0000_${paramsline}/${inpfile}_initial0000  run_0000_${paramsline}/particles_break_index_${breaktime} 0 > segment_restart_iter0.log
./move_files.sh segment_0000_${paramsline} ${inpfile}.inp
#
mkdir geometry
cp *.xyz geometry/
cp *.hdf5 geometry/
cp *initial*.xdmf geometry/
cp *initial*.h5 geometry/
#
##./move_files.sh flexion_30s_gm4.0_kn5e4_bn2e5_5e-5_v0.105_equal_radii
##flexion_150s_maxt3.0s_gm4.0_kn5e2_bn1e4_1e-4_v0.025_equal_radii
##XXX CURRENT FILENAME ->flexion_3.0s_gm4.0_kn5e2_bn1e4_1e-4_v0.025_equal_radii :
##time ./simple_breaking_restart_equilv1 $filename 10 > simple_breaking_restart_kn${kn}_bn${bn}_equil.log#./generate_cylinder_logs.sh simple_breaking_restart_kn${kn}_bn${bn}_equil.log $maxt 1.0_kn${kn}_bn${bn}
##./move_files.sh flexion_${maxt}s_gm1.0_kn${kn}_bn${bn}_equal_radii_run_0000
##./cut_particles brick_restart_geometry 10 0 flexion_300s_maxt3.0s_gm8.0_kn5e2_bn1e4_1e-4_v0.025_equal_radii/brick_restart_test_initial0000 flexion_300s_maxt3.0s_gm8.0_kn5e2_bn1e4_1e-4_v0.025_equal_radii/particles_break_index_300.000100
##./move_files.sh flexion_300s_maxt3.0s_gm8.0_kn5e2_bn1e4_1e-4_v0.025_equal_radii_particlecut_0000

#Iterate through restart indices
for ((i=1 ; i<=$1 ; i++))
do
	echo "Running iteration: ${i}"
	#Load segmentation from previous step
	time ./simple_breaking_restart_equilv2 ${inpfile} 10 ${i} "load" segment_000$((i-1))_${paramsline}/brick_restart_geometry_initial000${i} > simple_breaking_restart_${paramsline}_iter${i}.log
	./generate_cylinder_logs.sh simple_breaking_restart_${paramsline}_iter${i}.log ${paramsline}
	cp *.xyz geometry/
	cp *.hdf5 geometry/
	cp *.pov geometry/
	cp brick_geometry_000${i}.* geometry/
	#cp *initial*.xdmf geometry/
	#cp *initial*.h5 geometry/
	./move_files.sh run_000${i}_${paramsline}
	./cut_particles brick_restart_geometry 10 0 run_000${i}_${paramsline}/${inpfile}_initial000${i}  run_000${i}_${paramsline}/particles_break_index_${breaktime} ${i} > segment_restart_iter${i}.log
	cp *.xyz geometry/
	cp *.hdf5 geometry/
	#cp *initial*.xdmf geometry/
	#cp *initial*.h5 geometry/
	./move_files.sh segment_000${i}_${paramsline} ${inpfile}.inp
done
