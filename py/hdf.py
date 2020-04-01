# coding: UTF-8
import numpy as np
import h5py

im = 512
jm = 512
km = 512
#im = 512
#jm = 512
#km = 512
num = 9
#num = 3
#num = 1
#num = 17
#nstep = 20
#itime = 1
#ftime = 104
itime = 1
ftime = 1
timejump=1
jumpsave=4
dir = '/Volumes/MAEDA_HDD-1/'
#dir = '/glv0/maedarn/dbug/'

for i in range(itime,ftime +1,timejump): #最後は含まない
    #f90 = open(dir+'Allmv/Allmv'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'dbug-ok/Alldphi'+"%03.f"%(i)+'.DAT', 'rb')
    f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'tag/tag'+"%03.f"%(i)+'001.DAT', 'rb')
    ary = np.fromfile(f90, np.float32,count=im*jm*km*num) #count=im*jm*dim*nstep
    #print(ary)
    uvhpy = ary.reshape(num,im,jm,km, order='F')
    #print(uvhpy)

    '''
    slice1=uvhpy[0,:,:,:]
    slice2=uvhpy[1,:,:,:]
    slice3=uvhpy[2,:,:,:]
    slice4=uvhpy[3,:,:,:]
    slice5=uvhpy[4,:,:,:]
    slice6=uvhpy[5,:,:,:]
    slice7=uvhpy[6,:,:,:]
    slice8=uvhpy[7,:,:,:]
    slice9=uvhpy[8,:,:,:]
    '''

    slice1=uvhpy[0,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice2=uvhpy[1,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice3=uvhpy[2,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    
    #slice1=uvhpy[0,0:im:4,0:jm:4,0:km:4]
    #slice2=uvhpy[1,0:im:4,0:jm:4,0:km:4]
    #slice3=uvhpy[2,0:im:4,0:jm:4,0:km:4]
    
    slice4=uvhpy[3,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice5=uvhpy[4,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice6=uvhpy[5,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice7=uvhpy[6,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice8=uvhpy[7,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    slice9=uvhpy[8,0:im:jumpsave,0:jm:jumpsave,0:km:jumpsave]
    


    """
    slice10=uvhpy[9,:,:,:]
    slice11=uvhpy[10,:,:,:]
    slice12=uvhpy[11,:,:,:]
    slice13=uvhpy[12,:,:,:]
    slice14=uvhpy[13,:,:,:]
    slice15=uvhpy[14,:,:,:]
    slice16=uvhpy[15,:,:,:]
    slice17=uvhpy[16,:,:,:]
    """
    #slice18=uvhpy[17,:,:,:]
    #slicejmp1=uvhpy[0,0:5:2,0:5:2,0:5:2]
    #slicejmp2=uvhpy[1,0:5:2,0:5:2,0:5:2]
    #import h5py

    with h5py.File(dir+'tag/dphiAllHDF'+"%03.f"%(i)+'.h5', 'w') as f:
      #  f.create_group('PHI')
      #  f['PHI'].create_dataset('PHI', data=slice1)
        '''
        f.create_group('tag')
        f['tag'].create_dataset('tag', data=slice1)
        
        f.create_group('den')
        f['den'].create_dataset('den', data=slice2)
        
        f.create_group('RHO')
        f['RHO'].create_dataset('RHO', data=slice3)
        '''
        
        
        f.create_group('RHO')
        f['RHO'].create_dataset('RHO', data=slice1)
        
        #f['hoge'].create_dataset('a', data=slice2)

        #f.create_group('jump')
        #f['jump'].create_dataset('xj', data=slicejmp1)
        #f['jump'].create_dataset('aj', data=slicejmp2)
        
        

        f.create_group('V')
        f['V'].create_dataset('VX', data=slice2)
        f['V'].create_dataset('VY', data=slice3)
        f['V'].create_dataset('VZ', data=slice4)

        f.create_group('P')
        f['P'].create_dataset('P', data=slice5)

        f.create_group('B')
        f['B'].create_dataset('BX', data=slice6)
        f['B'].create_dataset('BY', data=slice7)
        f['B'].create_dataset('BZ', data=slice8)

        f.create_group('PHI')
        f['PHI'].create_dataset('PHI', data=slice9)
        
        """
        f.create_group('RHO-h')
        f['RHO-h'].create_dataset('RHO-h', data=slice9)

        f.create_group('RHO-p')
        f['RHO-p'].create_dataset('RHO-p', data=slice10)

        f.create_group('RHO-h2')
        f['RHO-h2'].create_dataset('RHO-h2', data=slice11)

        f.create_group('RHO-he')
        f['RHO-he'].create_dataset('RHO-he', data=slice12)

        f.create_group('RHO-hep')
        f['RHO-hep'].create_dataset('RHO-hep', data=slice13)

        f.create_group('RHO-c')
        f['RHO-c'].create_dataset('RHO-c', data=slice14)

        f.create_group('RHO-co')
        f['RHO-co'].create_dataset('RHO-co', data=slice15)

        f.create_group('RHO-cp')
        f['RHO-cp'].create_dataset('RHO-cp', data=slice16)

        f.create_group('PHI')
        f['PHI'].create_dataset('PHI', data=slice17)

        """
