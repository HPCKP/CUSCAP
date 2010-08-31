#!/bin/bash
#set -xv

mpirun -n 1 /home/jblasco/src/VASP/vasp.5.2/vasp >> vasp.log &

sleep 3
ompi-checkpoint --status $(pidof mpirun)

