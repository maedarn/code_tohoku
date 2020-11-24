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
timejump=4
Msun=0.0240
dxyz=100.0/512.0
mH=1.0
mHe=4.0
mH2=2.0
mC=12.0
mCO=28.0

#dir = '/glv0/maedarn/clst-form-HIcol/DTF/DTFm10-disp0-2-nolog/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wsb-wg-sm/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wsb-wg-sm-200pc/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/Clmpanly/'
#dir = '/glv0/maedarn/dbug/debuggausssidelwb/'
dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-optthick/'
#dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-hiden10/'

for i in range(itime,ftime +1,timejump): #最後は含まない
    #f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
    f90 = open(dir+'All/All'+"%03.f"%(i)+'.DAT', 'rb')
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
    sliceshock = np.where(sliceshock > 1, 1, 0)
    
    slice1shockrho10 = np.where(slice1 >= 12.7, 1, 0)
    slice1shockrho50 = np.where(slice1 >= 12.7*5.0, 1, 0)
    slice1shockrho100 = np.where(slice1 >= 127.0, 1, 0)
    slice1shockrho500 = np.where(slice1 >= 127.0*5.0, 1, 0)
    slice1shockrho1000 = np.where(slice1 >= 1270.0, 1, 0)
    slice1shockrho5000 = np.where(slice1 >= 1270.0*5.0, 1, 0)
    slice1shockrho10000 = np.where(slice1 >= 12700.0, 1, 0)
    
    slice1shockrho30 = np.where(slice1 >= 12.7*3.1622, 1, 0)
    slice1shockrho300 = np.where(slice1 >= 127.0*3.1622, 1, 0)
    slice1shockrho3000 = np.where(slice1 >= 1270.0*3.1622, 1, 0)

    

    slice1shock10=sliceshock*slice1*slice1shockrho10
    #slice2shock=sliceshock*slice2
    #slice3shock=sliceshock*slice3
    #slice4shock=sliceshock*slice4
    #slice5shock=sliceshock*slice5
    #slice6shock=sliceshock*slice6
    #slice7shock=sliceshock*slice7
    #slice8shock=sliceshock*slice8
    Hshock10=sliceshock*slice9*slice1shockrho10
    protonshock10=sliceshock*slice10*slice1shockrho10
    H2shock10=sliceshock*slice11*slice1shockrho10
    Heshock10=sliceshock*slice12*slice1shockrho10
    Hepshock10=sliceshock*slice13*slice1shockrho10
    Cshock10=sliceshock*slice14*slice1shockrho10
    COshock10=sliceshock*slice15*slice1shockrho10
    Cpshock10=sliceshock*slice16*slice1shockrho10
    
    slice1shock50=sliceshock*slice1*slice1shockrho50
    Hshock50=sliceshock*slice9*slice1shockrho50
    protonshock50=sliceshock*slice10*slice1shockrho50
    H2shock50=sliceshock*slice11*slice1shockrho50
    Heshock50=sliceshock*slice12*slice1shockrho50
    Hepshock50=sliceshock*slice13*slice1shockrho50
    Cshock50=sliceshock*slice14*slice1shockrho50
    COshock50=sliceshock*slice15*slice1shockrho50
    Cpshock50=sliceshock*slice16*slice1shockrho50
    
    Hshock100=sliceshock*slice9*slice1shockrho100
    protonshock100=sliceshock*slice10*slice1shockrho100
    H2shock100=sliceshock*slice11*slice1shockrho100
    Heshock100=sliceshock*slice12*slice1shockrho100
    Hepshock100=sliceshock*slice13*slice1shockrho100
    Cshock100=sliceshock*slice14*slice1shockrho100
    COshock100=sliceshock*slice15*slice1shockrho100
    Cpshock100=sliceshock*slice16*slice1shockrho100
    
    slice1shock500=sliceshock*slice1*slice1shockrho500
    Hshock500=sliceshock*slice9*slice1shockrho500
    protonshock500=sliceshock*slice10*slice1shockrho500
    H2shock500=sliceshock*slice11*slice1shockrho500
    Heshock500=sliceshock*slice12*slice1shockrho500
    Hepshock500=sliceshock*slice13*slice1shockrho500
    Cshock500=sliceshock*slice14*slice1shockrho500
    COshock500=sliceshock*slice15*slice1shockrho500
    Cpshock500=sliceshock*slice16*slice1shockrho500
    
    Hshock1000=sliceshock*slice9*slice1shockrho1000
    protonshock1000=sliceshock*slice10*slice1shockrho1000
    H2shock1000=sliceshock*slice11*slice1shockrho1000
    Heshock1000=sliceshock*slice12*slice1shockrho1000
    Hepshock1000=sliceshock*slice13*slice1shockrho1000
    Cshock1000=sliceshock*slice14*slice1shockrho1000
    COshock1000=sliceshock*slice15*slice1shockrho1000
    Cpshock1000=sliceshock*slice16*slice1shockrho1000
    
    slice1shock5000=sliceshock*slice1*slice1shockrho5000
    Hshock5000=sliceshock*slice9*slice1shockrho5000
    protonshock5000=sliceshock*slice10*slice1shockrho5000
    H2shock5000=sliceshock*slice11*slice1shockrho5000
    Heshock5000=sliceshock*slice12*slice1shockrho5000
    Hepshock5000=sliceshock*slice13*slice1shockrho5000
    Cshock5000=sliceshock*slice14*slice1shockrho5000
    COshock5000=sliceshock*slice15*slice1shockrho5000
    Cpshock5000=sliceshock*slice16*slice1shockrho5000
    
    Hshock10000=sliceshock*slice9*slice1shockrho10000
    protonshock10000=sliceshock*slice10*slice1shockrho10000
    H2shock10000=sliceshock*slice11*slice1shockrho10000
    Heshock10000=sliceshock*slice12*slice1shockrho10000
    Hepshock10000=sliceshock*slice13*slice1shockrho10000
    Cshock10000=sliceshock*slice14*slice1shockrho10000
    COshock10000=sliceshock*slice15*slice1shockrho10000
    Cpshock10000=sliceshock*slice16*slice1shockrho10000
    
    slice1shock30=sliceshock*slice1*slice1shockrho30
    Hshock30=sliceshock*slice9*slice1shockrho30
    protonshock30=sliceshock*slice10*slice1shockrho30
    H2shock30=sliceshock*slice11*slice1shockrho30
    Heshock30=sliceshock*slice12*slice1shockrho30
    Hepshock30=sliceshock*slice13*slice1shockrho30
    Cshock30=sliceshock*slice14*slice1shockrho30
    COshock30=sliceshock*slice15*slice1shockrho30
    Cpshock30=sliceshock*slice16*slice1shockrho30
    
    Hshock300=sliceshock*slice9*slice1shockrho300
    protonshock300=sliceshock*slice10*slice1shockrho300
    H2shock300=sliceshock*slice11*slice1shockrho300
    Heshock300=sliceshock*slice12*slice1shockrho300
    Hepshock300=sliceshock*slice13*slice1shockrho300
    Cshock300=sliceshock*slice14*slice1shockrho300
    COshock300=sliceshock*slice15*slice1shockrho300
    Cpshock300=sliceshock*slice16*slice1shockrho300
    
    Hshock3000=sliceshock*slice9*slice1shockrho3000
    protonshock3000=sliceshock*slice10*slice1shockrho3000
    H2shock3000=sliceshock*slice11*slice1shockrho3000
    Heshock3000=sliceshock*slice12*slice1shockrho3000
    Hepshock3000=sliceshock*slice13*slice1shockrho3000
    Cshock3000=sliceshock*slice14*slice1shockrho3000
    COshock3000=sliceshock*slice15*slice1shockrho3000
    Cpshock3000=sliceshock*slice16*slice1shockrho3000
    
    #Phishock=sliceshock*slice17
    #Tshock=sliceshock*slice18
    
    
    
    #ntot=np.sum(slice1shock)*dxyz*dxyz*dxyz
    ntotHI10=np.sum(Hshock10)*dxyz*dxyz*dxyz
    ntotH210=np.sum(H2shock10)*dxyz*dxyz*dxyz
    ntotp10=np.sum(protonshock10)*dxyz*dxyz*dxyz
    ntotCO10=np.sum(COshock10)*dxyz*dxyz*dxyz
    ntotCp10=np.sum(Cpshock10)*dxyz*dxyz*dxyz
    ntotC10=np.sum(Cshock10)*dxyz*dxyz*dxyz
    ntotHe10=np.sum(Heshock10+Hepshock10)*dxyz*dxyz*dxyz
    
    ntotHI100=np.sum(Hshock100)*dxyz*dxyz*dxyz
    ntotH2100=np.sum(H2shock100)*dxyz*dxyz*dxyz
    ntotp100=np.sum(protonshock100)*dxyz*dxyz*dxyz
    ntotCO100=np.sum(COshock100)*dxyz*dxyz*dxyz
    ntotCp100=np.sum(Cpshock100)*dxyz*dxyz*dxyz
    ntotC100=np.sum(Cshock100)*dxyz*dxyz*dxyz
    ntotHe100=np.sum(Heshock100+Hepshock100)*dxyz*dxyz*dxyz
    
    ntotHI1000=np.sum(Hshock1000)*dxyz*dxyz*dxyz
    ntotH21000=np.sum(H2shock1000)*dxyz*dxyz*dxyz
    ntotp1000=np.sum(protonshock1000)*dxyz*dxyz*dxyz
    ntotCO1000=np.sum(COshock1000)*dxyz*dxyz*dxyz
    ntotCp1000=np.sum(Cpshock1000)*dxyz*dxyz*dxyz
    ntotC1000=np.sum(Cshock1000)*dxyz*dxyz*dxyz
    ntotHe1000=np.sum(Heshock1000+Hepshock1000)*dxyz*dxyz*dxyz
    
    ntotHI10000=np.sum(Hshock10000)*dxyz*dxyz*dxyz
    ntotH210000=np.sum(H2shock10000)*dxyz*dxyz*dxyz
    ntotp10000=np.sum(protonshock10000)*dxyz*dxyz*dxyz
    ntotCO10000=np.sum(COshock10000)*dxyz*dxyz*dxyz
    ntotCp10000=np.sum(Cpshock10000)*dxyz*dxyz*dxyz
    ntotC10000=np.sum(Cshock10000)*dxyz*dxyz*dxyz
    ntotHe10000=np.sum(Heshock10000+Hepshock10000)*dxyz*dxyz*dxyz
    
    ntotHI50=np.sum(Hshock50)*dxyz*dxyz*dxyz
    ntotH250=np.sum(H2shock50)*dxyz*dxyz*dxyz
    ntotp50=np.sum(protonshock50)*dxyz*dxyz*dxyz
    ntotCO50=np.sum(COshock50)*dxyz*dxyz*dxyz
    ntotCp50=np.sum(Cpshock50)*dxyz*dxyz*dxyz
    ntotC50=np.sum(Cshock50)*dxyz*dxyz*dxyz
    ntotHe50=np.sum(Heshock50+Hepshock50)*dxyz*dxyz*dxyz
    
    ntotHI500=np.sum(Hshock500)*dxyz*dxyz*dxyz
    ntotH2500=np.sum(H2shock500)*dxyz*dxyz*dxyz
    ntotp500=np.sum(protonshock500)*dxyz*dxyz*dxyz
    ntotCO500=np.sum(COshock500)*dxyz*dxyz*dxyz
    ntotCp500=np.sum(Cpshock500)*dxyz*dxyz*dxyz
    ntotC500=np.sum(Cshock500)*dxyz*dxyz*dxyz
    ntotHe500=np.sum(Heshock500+Hepshock500)*dxyz*dxyz*dxyz
    
    ntotHI5000=np.sum(Hshock5000)*dxyz*dxyz*dxyz
    ntotH25000=np.sum(H2shock5000)*dxyz*dxyz*dxyz
    ntotp5000=np.sum(protonshock5000)*dxyz*dxyz*dxyz
    ntotCO5000=np.sum(COshock5000)*dxyz*dxyz*dxyz
    ntotCp5000=np.sum(Cpshock5000)*dxyz*dxyz*dxyz
    ntotC5000=np.sum(Cshock5000)*dxyz*dxyz*dxyz
    ntotHe5000=np.sum(Heshock5000+Hepshock5000)*dxyz*dxyz*dxyz
    
    ntotHI30=np.sum(Hshock30)*dxyz*dxyz*dxyz
    ntotH230=np.sum(H2shock30)*dxyz*dxyz*dxyz
    ntotp30=np.sum(protonshock30)*dxyz*dxyz*dxyz
    ntotCO30=np.sum(COshock30)*dxyz*dxyz*dxyz
    ntotCp30=np.sum(Cpshock30)*dxyz*dxyz*dxyz
    ntotC30=np.sum(Cshock30)*dxyz*dxyz*dxyz
    ntotHe30=np.sum(Heshock30+Hepshock30)*dxyz*dxyz*dxyz
    
    ntotHI300=np.sum(Hshock300)*dxyz*dxyz*dxyz
    ntotH2300=np.sum(H2shock300)*dxyz*dxyz*dxyz
    ntotp300=np.sum(protonshock300)*dxyz*dxyz*dxyz
    ntotCO300=np.sum(COshock300)*dxyz*dxyz*dxyz
    ntotCp300=np.sum(Cpshock300)*dxyz*dxyz*dxyz
    ntotC300=np.sum(Cshock300)*dxyz*dxyz*dxyz
    ntotHe300=np.sum(Heshock300+Hepshock300)*dxyz*dxyz*dxyz
    
    ntotHI3000=np.sum(Hshock3000)*dxyz*dxyz*dxyz
    ntotH23000=np.sum(H2shock3000)*dxyz*dxyz*dxyz
    ntotp3000=np.sum(protonshock3000)*dxyz*dxyz*dxyz
    ntotCO3000=np.sum(COshock3000)*dxyz*dxyz*dxyz
    ntotCp3000=np.sum(Cpshock3000)*dxyz*dxyz*dxyz
    ntotC3000=np.sum(Cshock3000)*dxyz*dxyz*dxyz
    ntotHe3000=np.sum(Heshock3000+Hepshock3000)*dxyz*dxyz*dxyz
    
    

    
    stime=i*1.0
    #stime=i*0.25+0.25
    
    smean = [[0.0 for ii in range(8)] for jj in range(10)]
    
    #print(np.mean(a, axis=0))
    #print(np.mean(a, axis=1))
    
    #smean = [stime,Mtot,MtotHI,MtotH2,MtotCO,MtotHe,Mtotp,MtotCO,MtotCp,MtotC,ntotH,ntotCsp,ntotHI,ntotH2,ntotCO,ntotHe,ntotp,ntotCO,ntotCp,ntotC]
    smean[0[:] = [10.0,ntotHI10,ntotH210,ntotp10,ntotCO10,ntotCp10,ntotC10,ntotHe10]
    smean[1][:] = [30.0,ntotHI30,ntotH230,ntotp30,ntotCO30,ntotCp30,ntotC30,ntotHe30]
    smean[2][:] = [50.0,ntotHI50,ntotH250,ntotp50,ntotCO50,ntotCp50,ntotC50,ntotHe50]
    smean[3][:] = [100.0,ntotHI100,ntotH2100,ntotp100,ntotCO100,ntotCp100,ntotC100,ntotHe100]
    smean[4][:] = [300.0,ntotHI300,ntotH2300,ntotp300,ntotCO300,ntotCp300,ntotC300,ntotHe300]
    smean[5][:] = [500.0,ntotHI500,ntotH2500,ntotp500,ntotCO500,ntotCp500,ntotC500,ntotHe500]
    smean[6][:] = [1000.0,ntotHI1000,ntotH21000,ntotp1000,ntotCO1000,ntotCp1000,ntotC1000,ntotHe1000]
    smean[7][:] = [3000.0,ntotHI3000,ntotH23000,ntotp3000,ntotCO3000,ntotCp3000,ntotC3000,ntotHe3000]
    smean[8][:] = [5000.0,ntotHI5000,ntotH25000,ntotp5000,ntotCO5000,ntotCp5000,ntotC5000,ntotHe5000]
    smean[9][:] = [10000.0,ntotHI10000,ntotH210000,ntotp10000,ntotCO10000,ntotCp10000,ntotC10000,ntotHe10000]
    #smean = [stime,Mtot100,Mtot100HI,Mtot100H2,Mtot100HI/Mtot100,Mtot100H2/Mtot100]
    #smean = [stime,Mtot,MtotHI,MtotH2,MtotCO,MtotHe,MtotHI/Mtot,MtotH2/Mtot,MtotCO/Mtot,MtotHe/Mtot]
    #smean = [stime,s1,s5]
    
    for n in range(0,10,1):
        with open(dir+'AllHDF/MassratioMol-nsp.csv', mode='a') as file:
            arr=smean[n][:]
            writer = csv.writer(file, lineterminator='\n')
            writer.writerow(arr)
    
  

