#!/bin/sh
#------ pjsub option --------#
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-debug
#PJM -L node=1
#PJM --mpi proc=48
#PJM -L elapse=1:00:00
#PJM -j
#------- Program execution -------#
mpiexec ./a.out 
