#!/bin/bash -x
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-small
#PJM -L node=2
#PJM --mpi proc=4
#PJM -L elapse=1:00:00
#PJM -j
#PJM -S
module avail
module list
export OMP_NUM_THREADS=24
mpiexec ./a.out
