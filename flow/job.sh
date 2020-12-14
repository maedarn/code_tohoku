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

#!/bin/sh
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-small
#PJM -L node=2
#PJM --mpi proc=64
#PJM -L elapse=1:00:00
#PJM -j
#PJM -S
#PJM --mail-list maeda.ryunosuke@j.mbox.nagoya-u.ac.jp
#PJM -m b
#PJM -m e
#PJM -o log.txt
#PJM -N grvwv-test-1

fipp -C -d./tmp mpiexec -n 64 ./a.out
