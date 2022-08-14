#!/bin/bash

grep "Upper cylinder" $1 > upper_cylinder_$2_gm$3_30s.log
grep "Lower cylinder a" $1 > lower_cylinder_a_$2_gm$3_30s.log
grep "Lower cylinder b" $1 > lower_cylinder_b_$2_gm$3_30s.log
grep "Upper particle" $1 > upper_particle_$2_gm$3_30s.log
