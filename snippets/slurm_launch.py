#!/usr/bin/env python3

""" launch a lot of jobs to study the performance of the initialiser """

def launchBestOfN ( bestOfN : int ):
    """ launch best-of-<n> """
    repetitions = int ( 1000. / bestOfN )
    basedir = f"/scratch-cbe/users/{os.environ['USER']}"
    outputdir = f"{basedir}/outputs"
    cmd = f"sbatch --time 479 --error {outputdir}/init{bestOfN}.out --output {outputdir}/init{bestOfN}.out ./crIS{bestOfN}.sh"
    for i in range(repetitions):
        o = subprocess.getoutput ( cmd )
        print ( cmd, o )

def launchAll():
    import subprocess
    seq = [ 3, 5, 10, 20, 50, 100, 200 ]
    for bestOfN in seq[::-1]:
        launchBestOfN ( bestOfN )

if __name__ == "__main__":
    launchAll()
