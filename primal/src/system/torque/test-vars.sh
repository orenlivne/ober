#!/bin/bash
# Test oassing env vars to a PBS script. Run with "qsub -v VAR1=1,VAR2=2 test-vars.sh"
echo $HOSTNAME
echo "VAR1=$VAR1"
echo "VAR2=$VAR2"
