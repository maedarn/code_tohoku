#!/bin/sh
#------ pjsub option --------#
#PJM -g z40309n 
#PJM -N Zm4V100D1B1A15
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-small
#PJM -L node=11
#PJM --mpi proc=512
#PJM -L elapse=5:00:00
#PJM -j
#PJM -o log.txt
#PJM -e Err.txt
#PJM -S
#PJM --mail-list maeda.ryunosuke@j.mbox.nagoya-u.ac.jp
#PJM -mber
#------- Program execution -------#

cd /data/group1/z40309n/lmt/Zm4V100D1B1A45/
mpiexec ./Zm4V100D1B1A45.out
