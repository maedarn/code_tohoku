#!/usr/bin/bash
#echo "test"が10回実行される。
#bashの変数は全てを文字列として扱う

user="maedarn"
host="xc.cfca.nao.ac.jp"
dicr="/work/maedarn/3DMHD/samplecnv2/"

gfortran kaiseki3D-loop.f90

for i in `seq -f %03g 142 186` #1-10の配列
do
    #echo "test"
    scp ${user}@${host}:${dicr}${i}*.dat .
#    gfortran kaiseki3.f90
    ./a.out
#    i
    rm ${i}*.dat
#    gfortran test.f90
#    ./a.out #<<
    # ${i}
    # EOS
#    echo $i
#    ${i} -> 'a.out'
done
