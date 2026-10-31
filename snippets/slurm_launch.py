#!/usr/bin/env python3

""" launch a lot of jobs to study the performance of the initialiser """
    
import os, subprocess

SQUEUE_FORMAT = "%.8i %.9P %.80j %.8u %.10T %.10M %.6D %R"
SQUEUE = f'squeue --me -o "{SQUEUE_FORMAT}" --sort=-i'

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
    cmd = f"for i in `{SQUEUE} | grep cr  | cut -d ' ' -f1`; do scancel $i; done"
    subprocess.getoutput ( cmd )


def launchAll():
    seq = [ 3, 5, 10, 20, 50, 100, 200 ]
    for bestOfN in seq[::-1]:
        launchBestOfN ( bestOfN )

if __name__ == "__main__":
    launchAll()
