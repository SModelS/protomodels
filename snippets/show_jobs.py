#!/usr/bin/env python3

import subprocess

def showOneJob( jobid : int ):
    cmd = f"jobinfo {jobid}" 
    output = subprocess.getoutput ( cmd ).split("\n")
    for line in output:
        if "State" in line:
            token = line.replace("State","")
            token = token.strip()
            token = token[2:]
            print ( f"{jobid} {token}" )

def show():
    with open ( "jobs", "rt" ) as f:
        lines = f.readlines()
    for line in lines:
        jobid = line.strip()
        showOneJob ( int ( jobid ) )

if __name__ == "__main__":
    show()
