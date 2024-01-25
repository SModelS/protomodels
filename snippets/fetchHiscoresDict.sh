#!/bin/sh

# fetch hiscores dictionary file

rm -r hiscore.hi

scp clip-login-1:rundir/hiscores.dict .
