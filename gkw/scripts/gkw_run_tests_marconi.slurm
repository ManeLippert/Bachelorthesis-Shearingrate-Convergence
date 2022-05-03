#! /bin/bash
## To run the GKW tests on marconi, submit this job script:
##
##    sbatch -A <youraccount> ~/path/to/this/script
##
##
## One can submit this job from any folder, the $PWD does not need to be a
## certain place.
#
## This does NOT request the maximum resources possible on the debugging queue:
#SBATCH --partition=skl_fua_dbg
#SBATCH --nodes 2 #maximum would be 4
#SBATCH --time=2:00:00
#SBATCH --job-name GKW-tests

$GKW_HOME/scripts/gkw_run_tests -m > $GKW_HOME/tests.out
