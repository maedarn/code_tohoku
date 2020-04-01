# coding: UTF-8
import numpy as np
import h5py
import csv


#im = 64
#jm = 64
#km = 64
im = 128
jm = 128
km = 128
#im = 512
#jm = 512
#km = 512
#num = 1
#num = 3
#num = 7
#num = 17
num = 18
nstep = 20
#itime = 1
#ftime = 104
itime = 1
ftime = 100
timejump=1
#dir = '/glv0/maedarn/clst-form-HIcol/DTF/DTFm10-disp0-2-nolog/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wsb-wg-sm-200pc/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/Clmpanly/'
#dir = '/glv0/maedarn/dbug/debuggausssidelwb/'
dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-optthick/'

for i in range(itime,ftime +1,timejump): #最後は含まない
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    f90 = open(dir+'AllHDF/AllHDF'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'dbug-ok/Alldphi'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'tag/tag'+"%03.f"%(i)+'003.DAT', 'rb')
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'Clmpanly/CN098001.DAT', 'rb')
    #f90 = open(dir+'AllD.DAT', 'rb')
    ary = np.fromfile(f90, np.float32,count=im*jm*km*num) #count=im*jm*dim*nstep
    #ary = np.fromfile(f90, np.int32,count=im*jm*km*num)
    #print(ary)
    uvhpy = ary.reshape(num,im,jm,km, order='F')
    #print(uvhpy)
    
    #slice2=uvhpy[1,:,:,:]
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

    sliceshockrho = np.where(slice1 >= 10.0, 1, 0)
    sliceshockvx = np.where(abs(slice2) <= 40.0, 1, 0)
    
    sliceshock = sliceshockrho + sliceshockvx
    sliceshock = np.where(sliceshock >= 1, 1, 0)
    


    slice1shock=sliceshock*slice1
    slice2shock=sliceshock*slice2
    slice3shock=sliceshock*slice3
    slice4shock=sliceshock*slice4
    slice5shock=sliceshock*slice5
    slice6shock=sliceshock*slice6
    slice7shock=sliceshock*slice7
    slice8shock=sliceshock*slice8
    slice9shock=sliceshock*slice9
    slice10shock=sliceshock*slice10
    slice11shock=sliceshock*slice11
    slice12shock=sliceshock*slice12
    slice13shock=sliceshock*slice13
    slice14shock=sliceshock*slice14
    slice15shock=sliceshock*slice15
    slice16shock=sliceshock*slice16
    slice17shock=sliceshock*slice17
    slice18shock=sliceshock*slice18
    
    #a[a < 5].mean()
    s1=slice1shock[slice1shock>0].mean()
    s2=slice2shock[slice2shock!=0].std()
    s3=slice3shock[slice3shock!=0].std()
    s4=slice4shock[slice4shock!=0].std()
    s5=slice5shock[slice5shock>0].mean()
    s6=slice6shock[slice6shock!=0].mean()
    s7=slice7shock[slice7shock!=0].mean()
    s8=slice8shock[slice8shock!=0].mean()
    s9=slice9shock[slice9shock>0].mean()
    s10=slice10shock[slice10shock>0].mean()
    s11=slice11shock[slice11shock>0].mean()
    s12=slice12shock[slice12shock>0].mean()
    s13=slice13shock[slice13shock>0].mean()
    s14=slice14shock[slice14shock>0].mean()
    s15=slice15shock[slice15shock>0].mean()
    s16=slice16shock[slice16shock>0].mean()
    s17=slice17shock[slice17shock!=0].mean()
    s18=slice18shock[slice18shock>0].mean()
    stime=i*0.25
    
    #print(np.mean(a, axis=0))
    #print(np.mean(a, axis=1))
    
    smean = [stime,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18]
    #smean = [stime,s1,s5]
    
    
    with open(dir+'AllHDF/psrmeantimeallval.csv', mode='a') as file:
#f.write(smean)
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(smean)
    
    
"""
#print("2次元配列で、15以上の場合は1、それ以外は0に配列の要素を置換")
#upNpArray = np.where(npArray >= 15, 1, 0)
#print(upNpArray)

#print("2次元配列で、15以上の場合は元の配列のままとし、それ以外は0に配列の要素を置換")
#upNpArray = np.where(npArray >= 15, npArray, 0)
#print(upNpArray)
    with h5py.File(dir+'AllHDF/PSRHDF'+"%03.f"%(i)+'.h5', 'w') as f:
    #  f.create_group('PHI')
    #  f['PHI'].create_dataset('PHI', data=slice1)
    
        f.create_group('RHO')
        f['RHO'].create_dataset('RHO', data=slice1)

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

        f.create_group('RHOpsr')
        f['RHOpsr'].create_dataset('RHOpsr', data=slice1shock)

        f.create_group('Vpsr')
        f['Vpsr'].create_dataset('VXpsr', data=slice2shock)
        f['Vpsr'].create_dataset('VYpsr', data=slice3shock)
        f['Vpsr'].create_dataset('VZpsr', data=slice4shock)

        f.create_group('Ppsr')
        f['Ppsr'].create_dataset('Ppsr', data=slice5shock)
    
        f.create_group('Bpsr')
        f['Bpsr'].create_dataset('BXpsr', data=slice6shock)
        f['Bpsr'].create_dataset('BYpsr', data=slice7shock)
        f['Bpsr'].create_dataset('BZpsr', data=slice8shock)
"""
