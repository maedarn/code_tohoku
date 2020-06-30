# coding: UTF-8
import numpy as np
import h5py


#im = 64
#jm = 64
#km = 64
im = 128
jm = 128
km = 128
#im = 512
#jm = 512
#km = 512
num = 25
nstep = 20
#itime = 1
#ftime = 104
itime = 0
ftime = 100
timejump=1
dir = '/glv0/maedarn/test-grvwave/PHIINI/'

for i in range(itime,ftime +1,timejump): #最後は含まない
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'dbug-ok/Alldphi'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'tagn/tag'+"%03.f"%(i)+'001.DAT', 'rb')
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'Clmpanly/CN098001.DAT', 'rb')
    #f90 = open(dir+'AllD.DAT', 'rb')
    ary = np.fromfile(f90, np.float32,count=im*jm*km*num) #count=im*jm*dim*nstep
    #ary = np.fromfile(f90, np.int32,count=im*jm*km*num)
    #print(ary)
    uvhpy = ary.reshape(num,im,jm,km, order='F')
    
    slice1=uvhpy[0,:,:,:]
    slice2=uvhpy[1,:,:,:]
    slice3=uvhpy[2,:,:,:]
    slice4=uvhpy[3,:,:,:]
    slice5=uvhpy[4,:,:,:]
    slice6=uvhpy[5,:,:,:]
    slice7=uvhpy[6,:,:,:]
    slice8=uvhpy[7,:,:,:]
    slice9=uvhpy[8,:,:,:]
    slice10=uvhpy[9,:,:,:]
    slice11=uvhpy[10,:,:,:]
    slice12=uvhpy[11,:,:,:]
    slice13=uvhpy[12,:,:,:]
    slice14=uvhpy[13,:,:,:]
    slice15=uvhpy[14,:,:,:]
    slice16=uvhpy[15,:,:,:]
    slice17=uvhpy[16,:,:,:]
    slice18=uvhpy[17,:,:,:]
    slice19=uvhpy[18,:,:,:]
    slice20=uvhpy[19,:,:,:]
    slice21=uvhpy[20,:,:,:]
    slice22=uvhpy[21,:,:,:]
    slice23=uvhpy[22,:,:,:]
    slice24=uvhpy[23,:,:,:]
    slice25=uvhpy[24,:,:,:]

    """
    slice1=uvhpy[0,0:im:4,0:jm:4,0:km:4]
    slice2=uvhpy[1,0:im:4,0:jm:4,0:km:4]
    slice3=uvhpy[2,0:im:4,0:jm:4,0:km:4]
    slice4=uvhpy[3,0:im:4,0:jm:4,0:km:4]
    """
    #slice5=uvhpy[4,:,:,:]    
    #slice5=uvhpy[4,0:im:4,0:jm:4,0:km:4]
    '''    
    slice6=uvhpy[5,0:im:4,0:jm:4,0:km:4]
    slice7=uvhpy[6,0:im:4,0:jm:4,0:km:4]
    '''
    '''
    slice8=uvhpy[7,0:im:4,0:jm:4,0:km:4]
    slice9=uvhpy[8,0:im:4,0:jm:4,0:km:4]
    '''
    '''
    slice10=uvhpy[9,0:im:4,0:jm:4,0:km:4]
    slice11=uvhpy[10,0:im:4,0:jm:4,0:km:4]
    slice12=uvhpy[11,0:im:4,0:jm:4,0:km:4]
    slice13=uvhpy[12,0:im:4,0:jm:4,0:km:4]
    slice14=uvhpy[13,0:im:4,0:jm:4,0:km:4]
    slice15=uvhpy[14,0:im:4,0:jm:4,0:km:4]
    slice16=uvhpy[15,0:im:4,0:jm:4,0:km:4]
    slice17=uvhpy[16,0:im:4,0:jm:4,0:km:4]
    
    slice18=uvhpy[17,0:im:4,0:jm:4,0:km:4]
    '''
    #slicejmp1=uvhpy[0,0:5:2,0:5:2,0:5:2]
    #slicejmp2=uvhpy[1,0:5:2,0:5:2,0:5:2]
    #import h5py

    #with h5py.File(dir+'AllHDF/HDF'+"%03.f"%(i)+'.h5', 'w') as f:
    with h5py.File(dir+'All/AllHDF'+"%03.f"%(i)+'.h5', 'w') as f:
        f.create_group('Phiwv')
        f['Phiwv'].create_dataset('Phiwv1', data=slice1)
        f['Phiwv'].create_dataset('Phiwv2', data=slice2)
        f['Phiwv'].create_dataset('Phiwv3', data=slice3)
        f['Phiwv'].create_dataset('Phiwv4', data=slice4)
        f['Phiwv'].create_dataset('Phiwv5', data=slice5)
        f['Phiwv'].create_dataset('Phiwv6', data=slice6)
        f['Phiwv'].create_dataset('Phiwv7', data=slice7)
        f['Phiwv'].create_dataset('Phiwv8', data=slice8)
        
        f.create_group('Phigrdwv')
        f['Phigrdwv'].create_dataset('Phigrdwv1', data=slice9)
        f['Phigrdwv'].create_dataset('Phigrdwv2', data=slice10)
        f['Phigrdwv'].create_dataset('Phigrdwv3', data=slice11)
        f['Phigrdwv'].create_dataset('Phigrdwv4', data=slice12)
        f['Phigrdwv'].create_dataset('Phigrdwv5', data=slice13)
        f['Phigrdwv'].create_dataset('Phigrdwv6', data=slice14)
        f['Phigrdwv'].create_dataset('Phigrdwv7', data=slice15)
        f['Phigrdwv'].create_dataset('Phigrdwv8', data=slice16)

        f.create_group('Phiexa')
        f['Phiexa'].create_dataset('Phiexa', data=slice17)
        
        f.create_group('Phigrdwvexa')
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa1', data=slice17)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa2', data=slice18)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa3', data=slice19)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa4', data=slice20)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa5', data=slice21)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa6', data=slice22)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa7', data=slice23)
        f['Phigrdwvexa'].create_dataset('Phigrdwvexa8', data=slice24)
        
