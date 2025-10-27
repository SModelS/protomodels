#!/usr/bin/env python3

import os

def reverse_readline(filename, buf_size=1024):
    """Generator to read a file line by line in reverse order."""
    with open(filename, 'rb') as f:
        f.seek(0, os.SEEK_END)
        buffer = b""
        pos = f.tell()
        while pos > 0:
            read_size = min(buf_size, pos)
            pos -= read_size
            f.seek(pos)
            chunk = f.read(read_size)
            buffer = chunk + buffer
            *lines, buffer = buffer.split(b"\n")
            for line in reversed(lines):
                yield line.decode("utf-8")
        if buffer:
            yield buffer.decode("utf-8")

def showWalkerid ( walkerid : int ) -> tuple:
    filename = f"logs/walker_{walkerid}.log"
    nfin, ntot = 0, 0
    for line in reverse_readline( filename ):
        p1 = line.find ( "Step ")
        p2 = line.find ( "finished" )
        if p1>-1 and p2>-1:
            token = line[p1+5:p2]
            nfin_i, ntot_i = tuple(map(int,token.split("/")))
            nfin +=  nfin_i
            ntot += ntot_i
            print ( f"#{walkerid:2d}: {nfin_i:5d}/{ntot_i:5d} finished" )
            break
    return nfin, ntot

def show():
    import glob
    walks = glob.glob ( "logs/walker_*.log" )
    walkerids = set()
    for walk in walks:
        walkerid = walk.replace("logs/walker_","").replace(".log","")
        try:
            walkerids.add ( int( walkerid) )
        except ValueError as e:
            pass
    nfin, ntot = [], []
    for walkerid in walkerids:
        nfin_i, ntot_i = showWalkerid ( walkerid )
        nfin.append ( nfin_i )
        ntot.append ( ntot_i )
    perc = sum(nfin)/sum(ntot)*100.
    print ( f"total: {sum(nfin)}/{sum(ntot)} ({perc:.1f}%)" )

if __name__ == "__main__":
    show()
