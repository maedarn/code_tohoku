#!/usr/bin/bash
#echo "test"が10回実行される。
#bashの変数は全てを文字列として扱う

user="maedarn"
host="an.cfca.nao.ac.jp"
dicr="/home/maedarn/cnv100bgmldrd/Allnew"
al="Allnew"

gfortran kaiseki3Dforbrto.f90 -o test.out
cc -lhdf5 hdfcorrectfinal2.c

for i in `seq -f %03g 142 184` #1-10の配列
do
#echo "test"
scp ${user}@${host}:${dicr}${i}*.DAT .
./test.out
./a.out
#    i
rm ${al}${i}*.DAT
rm f${i}*.dat
#    gfortran test.f90
#    ./a.out #<<
# ${i}
# EOS
#    echo $i
#    ${i} -> 'a.out'
done

