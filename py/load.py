import numpy as np

"""
from scipy.io import FortranFile

im = 10
jm = 10
nstep = 10

f90 = open('test.dat', 'rb')
ary = np.fromfile(f90, np.float32, count=nstep*im*jm*3)

uvhpy = ary.reshape(im, jm, 3, 10, order='F')  # byte×202*481*3=4×291486=1165944
print(uvhpy[:, :, 2, 0])

"""
#N=3  #1レコード番号あたりに格納されているデータの数。
#M=3  #レコードの総数。

"""
f=open('txt.dat','rb')
dty=np.dtype([('data','<'+str(N)+'f')])
chunk=np.fromfile(f,dtype=dty,count=M)
data=np.array([chunk[j]['data'] for j in range(0,M)])
#for j in range(0,M):
#    data[j,:]=chunk[j]['data']
"""

#input = np.fromfile(infile, '>f')
#print(data)


im = 3
jm = 3
#nstep = 20

# %% 0.fromfileを用いたファイル読み込み
"""
!多分 open(11, 云々 ,access='stream') でデータべた書きにしないと駄目なのでは？
!unformatted だけだとバイナリで書き出すがレコード区切りが入る。 stream は Fortran2003 以降の属性だが、昔風には direct access 形式でもまぁ普通は行けなくもない。
"""
f90 = open('txt.dat', 'rb')
ary = np.fromfile(f90, np.float32,count=im*jm) #count=im*jm*dim*nstep
print(ary)
# %% 1.読み込んだn番目の要素を取り出し、fortran形式の4次元配列に
uvhpy = ary.reshape(jm,im, order='F')  # 1要素あたりbyte×202*481*3=4×291486=165944byte,order='F':並び替えかた

print(uvhpy)
