#!/bin/bash
#PBS -N testarray
#PBS -l walltime=00:01:00
#PBS -l mppwidth=24
#PBS -S /bin/bash
#PBS -q development
#PBS -j oe
#PBS -V
#PBS -v VAR1=value1,VAR2=value2

aprun /lustre/beagle/ober/users/oren/ober/system/torque/test-array-jau-script.sh
