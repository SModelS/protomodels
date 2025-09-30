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

def showWalkerid ( walkerid ):
    filename = f"logs/walker_{walkerid}.log"
    for line in reverse_readline( filename ):
        p1 = line.find ( "Step ")
        p2 = line.find ( "finished" )
        if p1>-1 and p2>-1:
            print ( walkerid, line[p1+5:p2], "finished" )
            break
    return

def show():
    import glob
    walks = glob.glob ( "logs/walker_*.log" )
    walkerids = set()
    for walk in walks:
        walkerid = walk.replace("logs/walker_","").replace(".log","")
        walkerids.add ( int( walkerid) )
    for walkerid in walkerids:
        showWalkerid ( walkerid )

if __name__ == "__main__":
    show()
