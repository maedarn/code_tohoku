#!/bin/sh
#------ pjsub option --------#
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-debug
#PJM -L node=1
#PJM --mpi proc=8
#PJM -L elapse=1:00:00
#PJM -j
#PJM -S
#PJM -o log.txt
#PJM -N grv_tel_cnvrg
#------- Program execution -------#

#fapp -C -d./rep1 -Hevent=pa1 mpiexec -n 8 ./pfm.out
#fapp -C -d./rep2 -Hevent=pa2 mpiexec -n 8 ./pfm.out
#fapp -C -d./rep3 -Hevent=pa3 mpiexec -n 8 ./pfm.out
#fapp -C -d./rep4 -Hevent=pa4 mpiexec -n 8 ./pfm.out
#fapp -C -d./rep5 -Hevent=pa5 mpiexec -n 8 ./pfm.out

fipp -C -d./tmp1 mpiexec ./pfm.out
