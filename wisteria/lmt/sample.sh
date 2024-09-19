#!/bin/sh
#------ pjsub option --------#
#PJM -g wo24i054 
#PJM -L rscgrp=regular-o
#PJM -L node=11
#PJM --mpi proc=512
#PJM -L elapse=4:00:00
#PJM -N Zm4V100D1B1A15
#PJM --mail-list maeda.ryunosuke@j.mbox.nagoya-u.ac.jp
#PJM -mber
#PJM -j
#PJM -o log.txt
#PJM -e Err.txt
#------- Program execution -------#

cd /work/wo24i054/h43000/Zm4V100D1B1A15/
mpiexec -n 512  ./Zm4V100D1B1A15.out
