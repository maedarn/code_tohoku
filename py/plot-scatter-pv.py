# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py




dir = '/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-optthick/'
im=128
jm=128
km=128
num=18
nstep = 20
#itime = 1
#ftime = 104
itime = 98
ftime = 98
timejump=1
    
    
#fig, ax = plt.subplots()
fig = plt.figure()
#ax1 = fig.add_subplot(111)
ax = fig.add_subplot(1,1,1)
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
    
    #x1=slice1.reshape(-1,-1,)
    #y1=slice5.reshape(-1,-1,)
    x1=slice1/1.27
    y1=slice5*116.0
    

    #ax.set_xlabel('n')  # x軸ラベル
    #ax.set_ylabel('p/k_{B}')  # y軸ラベル
    #ax.set_title(r'$\sin(x)$ and $\cos(x)$') # グラフタイトル
    # ax.set_aspect('equal') # スケールを揃える
    #ax.grid()            # 罫線
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.plot(x1, y1)
    #matplotlib.pyplot.scatter
    plt.scatter(x1, y1, s=0.1)# norm=LogNorm())#, s=20)#, c=None, marker='o', cmap=None, norm=None,vmin=None, vmax=None, alpha=None, linewidths=None,verts=None, edgecolors=None, hold=None, data=None,**kwargs)
    ax.set_yscale('log')  # y軸をlogスケールで描く
    ax.set_xscale('log')  # x軸をlogスケールで描く
    ax.set_xlim([0.1, 1000000]) # x方向の描画範囲を指定
    ax.set_ylim([100, 10000000])    # y方向の描画範囲を指定
    #ax.plot(t, y1, color=c1, label=l1)
    #ax.plot(t, y2, color=c2, label=l2)
    #ax.plot(t, y3, color=c3, label=l3)
    #ax.plot(t, y4, color=c4, label=l4)
    #ax.legend(loc=0)    # 凡例
    #fig.tight_layout()  # レイアウトの設定
    # plt.savefig('hoge.png') # 画像の保存
    #plt.show()
                    
                    
    plt.savefig("foo.png")
                    
   
      
