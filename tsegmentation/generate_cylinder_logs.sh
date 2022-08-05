#!/bin/bash

grep "Upper cylinder" $1 > upper_cylinder_2.2_gm$2_30s.log
grep "Lower cylinder a" $1 > lower_cylinder_a_2.2_gm$2_30s.log
grep "Lower cylinder b" $1 > lower_cylinder_b_2.2_gm$2_30s.log
