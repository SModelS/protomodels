#!/bin/bash

for i in {0..255}; do
  printf "\e[48;5;${i}m%3d\e[0m " "$i"
  if (( (i + 1) % 16 == 0 )); then echo; fi
done
