#!/bin/bash

#make simple_breaking_restartv0
#make simple_breaking_write_points
#For a given number of times
#Run the simulation until it breaks
restart=0
#If restart is zero the ptype should be voronoi and restart should be 0
time ./simple_breaking_restartv0 brick_restart_base 6 > simple_breaking_restart_${restart}.log
#Move the files for visualization
bash move_files.sh run_000${restart}
cp initial_points_${restart}.xyz run_000${restart}
cp brick_restart_base.inp run_000${restart}/
#Look at the brick_planes_walls.res file for the highest strain Energy Field and use that timestep for breaking the geometry
#Put that number and change the reset number in the brick_restart_geometry.inp
#Generate new centers for the segmented particles
time ./simple_breaking_write_points brick_restart_geometry 1
bash move_files.sh segment_000${restart}
cp segment_000${restart}/brick_restart_geometry.* geometry/
cp brick_restart_geometry.inp segment_000${restart}/
#Change the reset number in brick_restart_base.inp and rerun the siulation
