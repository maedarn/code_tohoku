# coding: UTF-8
import numpy as np
import h5py

im = 140
num = 2

f90 = open('KIcoolbr.dat', 'rb')
ary = np.fromfile(f90, np.float32,count=im*num)

uvhpy = ary.reshape(num,im, order='F') 
slice1=uvhpy[0,:]
slice2=uvhpy[1,:]

with h5py.File('KIcool.h5', 'w') as f:
    f.create_group('RHO')
    f['RHO'].create_dataset('RHO', data=slice1)

    f.create_group('P')
    f['P'].create_dataset('P', data=slice2)
