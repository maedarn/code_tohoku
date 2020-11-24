# coding: UTF-8
import numpy as np
import h5py
import csv


#im = 64
#jm = 64
#km = 64
#im = 128
#jm = 128
#km = 128
im = 512
jm = 512
km = 512
#num = 1
#num = 3
#num = 7
#num = 17
num = 18
nstep = 20
#itime = 1
#ftime = 104
itime = 98
ftime = 98
timejump=1
G=0.00011142
L=100.0
dxyz=100.0/(512.0)
#dir = '/glv0/maedarn/clst-form-HIcol/DTF/DTFm10-disp0-2-nolog/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wsb-wg-sm-200pc/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/Clmpanly/'
#dir = '/glv0/maedarn/dbug/debuggausssidelwb/'
dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-optthick/'

for j in range(itime,ftime +1,timejump): #最後は含まない
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    f90 = open(dir+'All/All'+"%03.f"%(j)+'.DAT', 'rb')
    f91 = open(dir+'tag/tag'+"%03.f"%(j)+'001.DAT', 'rb')
    #f90 = open(dir+'dbug-ok/Alldphi'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'tag/tag'+"%03.f"%(i)+'003.DAT', 'rb')
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    #f90 = open(dir+'Clmpanly/CN098001.DAT', 'rb')
    #f90 = open(dir+'AllD.DAT', 'rb')
    ary = np.fromfile(f90, np.float32,count=im*jm*km*num) #count=im*jm*dim*nstep
    ary1 = np.fromfile(f91, np.int32,count=im*jm*km) #count=im*jm*dim*nstep
    #ary = np.fromfile(f90, np.int32,count=im*jm*km*num)
    #print(ary)
    uvhpy = ary.reshape(num,im,jm,km, order='F')
    tag = ary1.reshape(im,jm,km, order='F')
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
    
    
    '''
    with h5py.File(dir+'AllHDF/ASFHDF'+"%03.f"%(j)+'.h5', 'w') as f:
    #  f.create_group('PHI')
    #  f['PHI'].create_dataset('PHI', data=slice1)
        f.create_group('RHO')
        f['RHO'].create_dataset('RHO', data=slice1)
        #f['RHO'].create_dataset('RHO2', data=slice1)
         
        f.create_group('V') 
        f['V'].create_dataset('VX', data=slice2)
        f['V'].create_dataset('VY', data=slice3)
        f['V'].create_dataset('VZ', data=slice4)
    '''
    #shock trace
    sliceshockrho = np.where(slice1 >= 10.0, 1, 0)
    sliceshockvx = np.where(abs(slice2) <= 40.0, 1, 0)
    sliceshock = sliceshockrho + sliceshockvx
    sliceshock = np.where(sliceshock >= 1, 1, 0)
    
    #density trace
    slicerho10000 = np.where(slice1 >= 12700.0, 1, 0)
    slicerho1000 = np.where(slice1 >= 1270.0, 1, 0)
    slicerho100 = np.where(slice1 >= 127.0, 1, 0)
    
    #tag
    tagtrace=np.where(tag == 9, 1, 0)
    
    #sarface
    slice1shock=sliceshock*slice1
    slice1rho10000=sliceshock*slicerho10000
    slice1rho1000=sliceshock*slicerho1000
    slice1rho100=sliceshock*slicerho100
    sliceshocksf=sliceshock*slice1
    
    slice1rho10000xp=sliceshock*slicerho10000
    slice1rho1000xp=sliceshock*slicerho1000
    slice1rho100xp=sliceshock*slicerho100
    slice1rho10000xm=sliceshock*slicerho10000
    slice1rho1000xm=sliceshock*slicerho1000
    slice1rho100xm=sliceshock*slicerho100
    slice1rho10000yp=sliceshock*slicerho10000
    slice1rho1000yp=sliceshock*slicerho1000
    slice1rho100yp=sliceshock*slicerho100
    slice1rho10000ym=sliceshock*slicerho10000
    slice1rho1000ym=sliceshock*slicerho1000
    slice1rho100ym=sliceshock*slicerho100
    slice1rho10000zp=sliceshock*slicerho10000
    slice1rho1000zp=sliceshock*slicerho1000
    slice1rho100zp=sliceshock*slicerho100
    slice1rho10000zm=sliceshock*slicerho10000
    slice1rho1000zm=sliceshock*slicerho1000
    slice1rho100zm=sliceshock*slicerho100
    
    #u1=slice1rho100zm
    
    #with h5py.File(dir+'AllHDF/ASFHDF'+"%03.f"%(j)+'.h5', 'w') as f:
        
        #f.create_group('RHO')
        #f['RHO'].create_dataset('RHO', data=sliceshockrho.astype(np.float))
        #f['RHO'].create_dataset('RHO2', data=sliceshockrho.astype(np.float))
    
    for i in range(0,im):
        sliceshocksf[i,:,:]=sliceshock[i+1-(im)*int(1/(im-i)),:,:]+sliceshock[i-1-(im)*int(1/(im+i)),:,:]
        
        
        slice1rho10000xp[i,:,:] =slice1rho10000[i+1-(im)*int(1/(im-i)),:,:]
        slice1rho10000xm[i,:,:]=slice1rho10000[i-1+(im)*int(im/(im+i)),:,:]
        slice1rho10000yp[:,i,:] =slice1rho10000[:,i+1-(jm)*int(1/(jm-i)),:]
        slice1rho10000ym[:,i,:]=slice1rho10000[:,i-1+(jm)*int(jm/(jm+i)),:]
        slice1rho10000zp[:,:,i] =slice1rho10000[:,:,i+1-(km)*int(1/(km-i))]
        slice1rho10000zm[:,:,i]=slice1rho10000[:,:,i-1+(km)*int(km/(km+i))]
        
        
        slice1rho1000xp[i,:,:] =slice1rho1000[i+1-(im)*int(1/(im-i)),:,:]
        slice1rho1000xm[i,:,:]=slice1rho1000[i-1+(im)*int(im/(im+i)),:,:]
        slice1rho1000yp[:,i,:] =slice1rho1000[:,i+1-(jm)*int(1/(jm-i)),:]
        slice1rho1000ym[:,i,:]=slice1rho1000[:,i-1+(jm)*int(jm/(jm+i)),:]
        slice1rho1000zp[:,:,i] =slice1rho1000[:,:,i+1-(km)*int(1/(km-i))]
        slice1rho1000zm[:,:,i]=slice1rho1000[:,:,i-1+(km)*int(km/(km+i))]
        
        
        slice1rho100xp[i,:,:] =slice1rho100[i+1-(im)*int(1/(im-i)),:,:]
        slice1rho100xm[i,:,:]=slice1rho100[i-1+(im)*int(im/(im+i)),:,:]
        slice1rho100yp[:,i,:] =slice1rho100[:,i+1-(jm)*int(1/(jm-i)),:]
        slice1rho100ym[:,i,:]=slice1rho100[:,i-1+(jm)*int(jm/(jm+i)),:]
        slice1rho100zp[:,:,i] =slice1rho100[:,:,i+1-(km)*int(1/(km-i))]
        slice1rho100zm[:,:,i]=slice1rho100[:,:,i-1+(km)*int(km/(km+i))]
        
        

    sliceshocksf= np.where(sliceshocksf == 1, 1, 0)
    
    slice1rho10000sf=slice1rho10000xp +slice1rho10000xm+slice1rho10000yp +slice1rho10000ym+slice1rho10000zp +slice1rho10000zm+slice1rho10000
    slice1rho10000sf=np.where(slice1rho10000sf <= 6, 1, 0)
    #slice1rho10000sf=slice1rho10000sf*slice1rho10000*tag
    slice1rho10000sf=slice1rho10000sf*slice1rho10000
    
    slice1rho1000sf=slice1rho1000xp +slice1rho1000xm+slice1rho1000yp +slice1rho1000ym+slice1rho1000zp +slice1rho1000zm+slice1rho1000
    slice1rho1000sf=np.where(slice1rho1000sf <= 6, 1, 0)
    slice1rho1000sf=slice1rho1000sf*slice1rho1000
    
    slice1rho100sf=slice1rho100xp +slice1rho100xm+slice1rho100yp +slice1rho100ym+slice1rho100zp +slice1rho100zm+slice1rho100
    slice1rho100sf=np.where(slice1rho100sf <= 6, 1, 0)
    slice1rho100sf=slice1rho100sf*slice1rho100
        
    
    #Gravity
    loop10000=np.sum(slice1rho10000sf)
    #loop10000=int(np.sum(slice1rho10000sf))
    ii=np.full(loop10000,0)
    jj=np.full(loop10000,0)
    kk=np.full(loop10000,0)
    #ii=np.zeros(int(loop10000),dtype=np.int32)
    #jj=np.zeros(loop10000,dtype=np.int32)
    #kk=np.zeros(loop10000ldtype=np.int32)
    #ii=ii.astype(int)
    #jj=jj.astype(int)
    #kk=kk.astype(int)
    #kk.dtype = 'int32'
    numloop=0
    for nn in range(0,km):
        for mm in range(0,jm):
            for ll in range(0,im):
                a=slice1rho10000sf[ll,mm,nn]
                if a > 0:
                    ii[numloop]=ll
                    jj[numloop]=mm
                    kk[numloop]=nn
                    numloop=numloop+1
             
    iimax=np.max(ii)
    iimin=np.min(ii)
    jjmax=np.max(jj)
    jjmin=np.min(jj)
    kkmax=np.max(kk)
    kkmin=np.min(kk)
    
    print(iimax)
    print(iimin)
    
    sfgrv10000=np.zeros(loop10000)
    for gg in range(0,loop10000):
        moji=gg
        print(moji)
        #sfgrv10000[gg]=G*tag[ii[gg],jj[gg],kk[gg]]*slice1[ii[gg],jj[gg],kk[gg]]*slice1rho10000[ii[gg],jj[gg],kk[gg]]*dxyz*dxyz*dxyz/np.sqrt(0.00000000000000001)
        sfgrv10000[gg]=G*slice1[ii[gg],jj[gg],kk[gg]]*slice1rho10000[ii[gg],jj[gg],kk[gg]]*dxyz*dxyz*dxyz/np.sqrt(0.00000000000000001)
        for nn in range(iimin,iimax):
            for mm in range(jjmin,jjmax):
                for ll in range(kkmin,kkmax):
                    #if ll != ii[gg] and mm != jj[gg] and nn != kk[gg]:
                    #    sfgrv10000[gg]=sfgrv10000[gg]-G*slice1[ll,mm,nn]*slice1rho10000[ll,mm,nn]*dxyz*dxyz*dxyz/np.sqrt(float((ll-ii[gg])**2+(mm-jj[gg])**2+(nn-kk[gg])**2)*dxyz**2+0.00000000000000001)
                    #sfgrv10000[gg]=sfgrv10000[gg]-G*tag[ii[gg],jj[gg],kk[gg]]*slice1[ll,mm,nn]*slice1rho10000[ll,mm,nn]*dxyz*dxyz*dxyz/np.sqrt(float((ll-ii[gg])**2+(mm-jj[gg])**2+(nn-kk[gg])**2)*dxyz**2+0.00000000000000001)
                    sfgrv10000[gg]=sfgrv10000[gg]-G*slice1[ll,mm,nn]*slice1rho10000[ll,mm,nn]*dxyz*dxyz*dxyz/np.sqrt(float((ll-ii[gg])**2+(mm-jj[gg])**2+(nn-kk[gg])**2)*dxyz**2+0.00000000000000001)
                    
    
    Vesc=np.sqrt(-2.0*sfgrv10000)
    Vescmaen=Vesc.mean()
    print(Vescmaen)
    
    plt.scatter(x1, y1, s=0.1)
    ax.set_xlim([0.1, 1000000]) # x方向の描画範囲を指定
    ax.set_ylim([100, 10000000])    # y方向の描画範囲を指定
    plt.savefig("scatter.png")
    
    #a[a < 5].mean()
    #s1=slice1shock[slice1shock>0].mean()
    #s2=slice2shock[slice2shock!=0].std()
   
    
    #print(np.mean(a, axis=0))
    #print(np.mean(a, axis=1))
    
    #smean = [stime,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18]
    #smean = [stime,s1,s5]
    
    
    #with open(dir+'AllHDF/psrmeantimeallval.csv', mode='a') as file:
    #f.write(smean)
      #  writer = csv.writer(file, lineterminator='\n')
       # writer.writerow(smean)
    
    

#print("2次元配列で、15以上の場合は1、それ以外は0に配列の要素を置換")
#upNpArray = np.where(npArray >= 15, 1, 0)
#print(upNpArray)

#print("2次元配列で、15以上の場合は元の配列のままとし、それ以外は0に配列の要素を置換")
#upNpArray = np.where(npArray >= 15, npArray, 0)
#print(upNpArray)
    #with h5py.File(dir+'AllHDF/SFHDF'+"%03.f"%(j)+'.h5', 'w') as f:
    #  f.create_group('PHI')
    #  f['PHI'].create_dataset('PHI', data=slice1)
        #f.create_group('RHO')
        #f['RHO'].create_dataset('RHO', data=slice1)
        #f['RHO'].create_dataset('RHO2', data=slice1)
         
        #f.create_group('V') 
        #f['V'].create_dataset('VX', data=slice2)
        #f['V'].create_dataset('VY', data=slice3)
        #f['V'].create_dataset('VZ', data=slice4)
        #f.create_group('RHO')
        #f['RHO'].create_dataset('RHO', data=sliceshockrho.astype(np.float))
        #f['RHO'].create_dataset('RHO1', data=sliceshockrho.astype(np.float))
        #f.create_group('SF')
        #f['SF'].create_dataset('SK1', data=sliceshockrho.astype(np.float))
        #f['SF'].create_dataset('SK2', data=slice1rho100zm)
        #f['SF'].create_dataset('SK3', data=slice1rho100zm.astype(np.float))
        #f['SF'].create_dataset('SK', data=sliceshocksf.astype(np.float))
        #f['SF'].create_dataset('10000', data=slice1rho10000sf.astype(np.float))
        #f['SF'].create_dataset('1000', data=slice1rho1000sf.astype(np.float))
        #f['SF'].create_dataset('100', data=slice1rho100sf.astype(np.float))
        


