#!/usr/bin/env python3

import subprocess
from clip.cliphelpers import getJobStatus, readJobIds

def writeJobIds( jobids : set ):
    """ write the job ids from jobs files """
    from base.locker import lock, unlock
    lock ( "jobs" )
    with open ( "jobs", "wt" ) as f:
        for jobid in jobids:
            f.write ( f"{jobid}\n" )
        f.close()
    unlock ( "jobs" )

def showOneJob( jobid : int ):
    print ( f"{jobid} {getJobStatus(jobid)[jobid]}" )

def show():
    jobids = readJobIds()
    statuses = getJobStatus ( jobids )
    cleaned_jobs = set()
    for jobid, status in  statuses.items():
        print ( f"{jobid}: {status}" )
        if status not in [ "completed" ]:
            cleaned_jobs.add ( jobid )
    writeJobIds ( cleaned_jobs )

if __name__ == "__main__":
    show()
