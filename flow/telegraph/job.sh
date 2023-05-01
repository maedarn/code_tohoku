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

mpiexec -n 8 ./pfm.out
