#!/bin/bash
# rsync test
remote="olivne@bios.cri.uchicago.edu"
echo "test" > test
rsync -auvPn test ${remote}:/
rm -f test
ssh ${remote} "rm -f test"
