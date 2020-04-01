# coding: UTF-8
"""
    Pythonにおける文字列はasciiコードでエンコードされるため、日本語のようなマルチバイト文字があるとエラーとなってしまいます。
    日本語をプログラム中で扱う場合には、プログラムの先頭で文字コードを指定し、指定した文字コードでファイルを保存する必要があります。
    Pythonプログラムの文字コードを指定し日本語（UTF-8）を扱えるようにするには、プログラムの1行目もしくは2行目に次のように記述します。
    # coding: エンコーディング名
    UTF-8の他にも代表的な日本語のエンコードとしてShift_JIS、EUC-JIS、 ISO-2022-JPなどがあります。
"""
import numpy as np

im = 4
jm = 4
km = 4
num = 2
#nstep = 20

# %% 0.fromfileを用いたファイル読み込み
"""
!多分 open(11, 云々 ,access='stream') でデータべた書きにしないと駄目なのでは？
!unformatted だけだとバイナリで書き出すがレコード区切りが入る。 stream は Fortran2003 以降の属性だが、昔風には direct access 形式でもまぁ普通は行けなくもない。
"""
f90 = open('txt.dat', 'rb')
ary = np.fromfile(f90, np.float32,count=im*jm*km*num) #count=im*jm*dim*nstep
print(ary)
# %% 1.読み込んだn番目の要素を取り出し、fortran形式の4次元配列に
#uvhpy = ary.reshape(num,km,jm,im, order='F')  # 1要素あたりbyte×202*481*3=4×291486=165944byte,order='F':並び替えかた
uvhpy = ary.reshape(num,im,jm,km, order='F')
print(uvhpy)



import h5py

# a[1:8:2]          # インデックス1〜7の要素を2ステップで

slice1=uvhpy[0,:,:,:]
slice2=uvhpy[1,:,:,:]
slicejmp1=uvhpy[0,0:5:2,0:5:2,0:5:2]
slicejmp2=uvhpy[1,0:5:2,0:5:2,0:5:2]
with h5py.File('output_group.h5', 'w') as f:
    f.create_group('hoge')
    f['hoge'].create_dataset('x', data=slice1)
    f['hoge'].create_dataset('a', data=slice2)

    f.create_group('jump')
    f['jump'].create_dataset('xj', data=slicejmp1)
    f['jump'].create_dataset('aj', data=slicejmp2)
