#!/usr/bin/env python3

""" launch a lot of jobs to study the performance of the initialiser """
    
import os, subprocess

def launchBestOfN ( bestOfN : int ):
    """ launch best-of-<n> """
    repetitions = int ( bestOfN ) # best-of-200 gets 200 processes, etc
    basedir = f"/scratch-cbe/users/{os.environ['USER']}"
    outputdir = f"{basedir}/outputs"
    cmd = f"sbatch --time 479 --error {outputdir}/init{bestOfN}.out --output {outputdir}/init{bestOfN}.out ./crIS{bestOfN}.sh"
    for i in range(repetitions):
        o = subprocess.getoutput ( cmd )
        print ( cmd, o )

def cancel():
    cmd = "for i in `slurm q | grep cr  | cut -d ' ' -f1`; do scancel $i; done"
    subprocess.getoutput ( cmd )


def launchAll():
    seq = [ 3, 5, 10, 20, 50, 100, 200 ]
    for bestOfN in seq[::-1]:
        launchBestOfN ( bestOfN )

if __name__ == "__main__":
    launchAll()
