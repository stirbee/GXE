#!/bin/bash
if [  $# -eq 0 ]
  then
    name=par
  else
    name=$1
fi

time ./gmatrix < $name.dat > $name.G-matrix.lst

