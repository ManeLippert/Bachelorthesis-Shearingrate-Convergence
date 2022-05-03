#! /bin/bash
## To run the GKW tests on the ITER cluster, submit this job script:
##
##    sbatch ~/path/to/this/script
##
##
## One can submit this job from any folder, the $PWD does not need to be a
## certain place.
#
#SBATCH --nodes 2
#SBATCH --time=2:00:00
#SBATCH --job-name GKW-tests
#SBATCH --exclusive
#SBATCH --error myjob.err
#SBATCH -p gen10_ib

$GKW_HOME/scripts/gkw_run_tests -m > $GKW_HOME/tests.out
