#!/usr/bin/env python3

import subprocess
from clip.cliphelpers import getJobStatus, readJobIds

def showOneJob( jobid : int ):
    print ( f"{jobid} {getJobStatus(jobid)[jobid]}" )

def show():
    jobids = readJobIds()
    statuses = getJobStatus ( jobids )
    for jobid, status in  statuses.items():
        print ( f"{jobid}: {status}" )

if __name__ == "__main__":
    show()
