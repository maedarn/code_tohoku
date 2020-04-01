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

    #sliceshockrho = np.where(slice1 >= 10.0, 1, 0)
    #sliceshockvx = np.where(abs(slice2) <= 40.0, 1, 0)
    
    #sliceshock = sliceshockrho + sliceshockvx
    #sliceshock = np.where(sliceshock >= 1, 1, 0)
    


    #slice1shock=sliceshock*slice1
    #slice2shock=sliceshock*slice2
    #slice3shock=sliceshock*slice3
    #slice4shock=sliceshock*slice4
    #slice5shock=sliceshock*slice5
    #slice6shock=sliceshock*slice6
    #slice7shock=sliceshock*slice7
    #slice8shock=sliceshock*slice8
    #slice9shock=sliceshock*slice9
    #slice10shock=sliceshock*slice10
    #slice11shock=sliceshock*slice11
    #slice12shock=sliceshock*slice12
    #slice13shock=sliceshock*slice13
    #slice14shock=sliceshock*slice14
    #slice15shock=sliceshock*slice15
    #slice16shock=sliceshock*slice16
    #slice17shock=sliceshock*slice17
    #slice18shock=sliceshock*slice18
    
    #a[a < 5].mean()
    s1=slice1.mean(axis=2)
    s2=slice2.mean(axis=2)
    s3=slice3.mean(axis=2)
    s4=slice4.mean(axis=2)
    s5=slice5.mean(axis=2)
    s6=slice6.mean(axis=2)
    s7=slice7.mean(axis=2)
    s8=slice8.mean(axis=2)
    s9=slice9.mean(axis=2)
    s10=slice10.mean(axis=2)
    s11=slice11.mean(axis=2)
    s12=slice12.mean(axis=2)
    s13=slice13.mean(axis=2)
    s14=slice14.mean(axis=2)
    s15=slice15.mean(axis=2)
    s16=slice16.mean(axis=2)
    s17=slice17.mean(axis=2)
    s18=slice18.mean(axis=2)
    
    s1d1=s1.mean(axis=1)
    s1d2=s2.mean(axis=1)
    s1d3=s3.mean(axis=1)
    s1d4=s4.mean(axis=1)
    s1d5=s5.mean(axis=1)
    s1d6=s6.mean(axis=1)
    s1d7=s7.mean(axis=1)
    s1d8=s8.mean(axis=1)
    s1d9=s9.mean(axis=1)
    s1d10=s10.mean(axis=1)
    s1d11=s11.mean(axis=1)
    s1d12=s12.mean(axis=1)
    s1d13=s13.mean(axis=1)
    s1d14=s14.mean(axis=1)
    s1d15=s15.mean(axis=1)
    s1d16=s16.mean(axis=1)
    s1d17=s17.mean(axis=1)
    s1d18=s18.mean(axis=1)
    
    '''
    s1d1=s1d1.T
    s1d2=s1d2.T
    s1d3=s1d3.T
    s1d4=s1d4.T
    s1d5=s1d5.T
    s1d6=s1d6.T
    s1d7=s1d7.T
    s1d8=s1d8.T
    s1d9=s1d9.T
    s1d10=s1d10.T
    s1d11=s1d11.T
    s1d12=s1d12.T
    s1d13=s1d13.T
    s1d14=s1d14.T
    s1d15=s1d15.T
    s1d16=s1d16.T
    s1d17=s1d17.T
    s1d18=s1d18.T

    #stime=i*0.25
    
    #print(np.mean(a, axis=0))
    #print(np.mean(a, axis=1))
    
    a = np.arange(im).T
    '''
    #smean = [a,s1d1,s1d2,s1d3,s1d4,s1d5,s1d6,s1d7,s1d8,s1d9,s1d10,s1d11,s1d12,s1d13,s1d14,s1d15,s1d16,s1d17,s1d18]
    #smean = [stime,s1,s5]
    '''
    s1d1=s1d1
    s1d2=s1d2
    s1d3=s1d3
    s1d4=s1d4
    s1d5=s1d5
    s1d6=s1d6
    s1d7=s1d7
    s1d8=s1d8
    s1d9=s1d9
    s1d10=s1d10
    s1d11=s1d11
    s1d12=s1d12
    s1d13=s1d13
    s1d14=s1d14
    s1d15=s1d15
    s1d16=s1d16
    s1d17=s1d17
    s1d18=s1d18
    '''
    
    #stime=i*0.25
    
    #print(np.mean(a, axis=0))
    #print(np.mean(a, axis=1))
    
    #a = np.arange(im)
    #smean = [a,s1d1,s1d2,s1d3,s1d4,s1d5,s1d6,s1d7,s1d8,s1d9,s1d10,s1d11,s1d12,s1d13,s1d14,s1d15,s1d16,s1d17,s1d18]
    #smeant = smean.T
    
    np.savetxt(dir+'AllHDF/meanxdrho'+"%03.f"%(i)+'.txt', s1d1)
    np.savetxt(dir+'AllHDF/meanxdvx'+"%03.f"%(i)+'.txt', s1d2)
    np.savetxt(dir+'AllHDF/meanxdvy'+"%03.f"%(i)+'.txt', s1d3)
    np.savetxt(dir+'AllHDF/meanxdvz'+"%03.f"%(i)+'.txt', s1d4)
    np.savetxt(dir+'AllHDF/meanxdp'+"%03.f"%(i)+'.txt', s1d5)
    np.savetxt(dir+'AllHDF/meanxdbx'+"%03.f"%(i)+'.txt', s1d6)
    np.savetxt(dir+'AllHDF/meanxdby'+"%03.f"%(i)+'.txt', s1d7)
    np.savetxt(dir+'AllHDF/meanxdbz'+"%03.f"%(i)+'.txt', s1d8)
    np.savetxt(dir+'AllHDF/meanxdh'+"%03.f"%(i)+'.txt', s1d9)
    np.savetxt(dir+'AllHDF/meanxdhp'+"%03.f"%(i)+'.txt', s1d10)
    np.savetxt(dir+'AllHDF/meanxdh2'+"%03.f"%(i)+'.txt', s1d11)
    np.savetxt(dir+'AllHDF/meanxdHe'+"%03.f"%(i)+'.txt', s1d12)
    np.savetxt(dir+'AllHDF/meanxdHep'+"%03.f"%(i)+'.txt', s1d13)
    np.savetxt(dir+'AllHDF/meanxdC'+"%03.f"%(i)+'.txt', s1d14)
    np.savetxt(dir+'AllHDF/meanxdCO'+"%03.f"%(i)+'.txt', s1d15)
    np.savetxt(dir+'AllHDF/meanxdCp'+"%03.f"%(i)+'.txt', s1d16)
    np.savetxt(dir+'AllHDF/meanxdphi'+"%03.f"%(i)+'.txt', s1d17)
    np.savetxt(dir+'AllHDF/meanxdT'+"%03.f"%(i)+'.txt', s1d18)
    # with open(dir+'AllHDF/meanxd'+"%03.f"%(i)+'.csv', 'w') as file:
#f.write(smean)
#for p in range(1,im +1,1):
#      writer = csv.writer(file, lineterminator='\n')
#       writer.writerow(a[p],s1d1[p])
    
    
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
