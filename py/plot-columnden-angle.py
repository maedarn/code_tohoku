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
ax1 = fig.add_subplot(111)
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
    
    s1=slice1.sum(axis=0)
    s1=s1/1.27*0.78*3*10**18
    s9=slice9.sum(axis=0)
    s9=s9*0.78*3*10**18
    s11=slice11.sum(axis=0)
    s11=s11*0.78*3*10**18
    """
    s2=slice2.sum(axis=3)
    s3=slice3.sum(axis=3)
    s4=slice4.sum(axis=3)
    s5=slice5.sum(axis=3)
    s6=slice6.sum(axis=3)
    s7=slice7.sum(axis=3)
    s8=slice8.sum(axis=3)
    s9=slice9.sum(axis=3)
    s10=slice10.sum(axis=3)
    s11=slice11.sum(axis=3)
    s12=slice12.sum(axis=3)
    s13=slice13.sum(axis=3)
    s14=slice14.sum(axis=3)
    s15=slice15.sum(axis=3)
    s16=slice16.sum(axis=3)
    s17=slice17.sum(axis=3)
    s18=slice18.sum(axis=3)
    """
    
    bx = np.arange(im)
    by = np.arange(jm)

    #P=P*116.0
    #Rho=Rho/1.27

    #c1,c2,c3,c4 = "blue","green","red","black"      # 各プロットの色
    #l1,l2,l3,l4 = "sin","cos","abs(sin)","sin**2"   # 各ラベル

    x = np.arange(0, im) #x軸の描画範囲の生成。0から10まで0.05刻み。
    y = np.arange(0, jm) #y軸の描画範囲の生成。0から10まで0.05刻み。
    #ax.set_xlabel('n')  # x軸ラベル
    #ax.set_ylabel('p/k_{B}')  # y軸ラベル
    #ax.set_title(r'$\sin(x)$ and $\cos(x)$') # グラフタイトル
    # ax.set_aspect('equal') # スケールを揃える
    #ax.grid()            # 罫線
    #plt.yscale('log')
    #plt.xscale('log')
    #ax.set_xlim([0.01, 1000]) # x方向の描画範囲を指定
    #ax.set_ylim([1, 100000])    # y方向の描画範囲を指定
    #ax.plot(t, y1, color=c1, label=l1)
    #ax.plot(t, y2, color=c2, label=l2)
    #ax.plot(t, y3, color=c3, label=l3)
    #ax.plot(t, y4, color=c4, label=l4)
    #ax.legend(loc=0)    # 凡例
    #fig.tight_layout()  # レイアウトの設定
    # plt.savefig('hoge.png') # 画像の保存
    #plt.show()

    mappable0 = ax1.pcolormesh(x,y,s1, cmap='coolwarm', norm=LogNorm()) # ここがポイント！
    pp = fig.colorbar(mappable0, ax=ax1, orientation="vertical")
    #plt.pcolormesh(ax, ay, s1, cmap='hsv',norm=LogNorm)
    #pp=plt.colorbar (orientation="vertical")
    #pp.set_label("Label", fontname="Arial", fontsize=24) #カラーバーのラベル
    #pp.set_aspect('equal')
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig("columntot.png")
    
    mappable0 = ax1.pcolormesh(x,y,s9, cmap='coolwarm', norm=LogNorm()) # ここがポイント！
    pp = fig.colorbar(mappable0, ax=ax1, orientation="vertical")
    #plt.pcolormesh(ax, ay, s1, cmap='hsv',norm=LogNorm)
    #pp=plt.colorbar (orientation="vertical")
    #pp.set_label("Label", fontname="Arial", fontsize=24) #カラーバーのラベル
    #pp.set_aspect('equal')
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig("columnh1.png")
    
    mappable0 = ax1.pcolormesh(x,y,s11, cmap='coolwarm', norm=LogNorm()) # ここがポイント！
    pp = fig.colorbar(mappable0, ax=ax1, orientation="vertical")
    #plt.pcolormesh(ax, ay, s1, cmap='hsv',norm=LogNorm)
    #pp=plt.colorbar (orientation="vertical")
    #pp.set_label("Label", fontname="Arial", fontsize=24) #カラーバーのラベル
    #pp.set_aspect('equal')
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig("columnh2.png")
    #plt.xlabel('X', fontsize=24)
    #plt.ylabel('Y', fontsize=24)

    #    matplotlib.pyplot.scatter(Rho, P, s=20, c=None, marker='o', cmap=None, norm=None,
    #vmin=None, vmax=None, alpha=None, linewidths=None,
    #                         verts=None, edgecolors=None, hold=None, data=None,
    #                         **kwargs)

                              #    matplotlib.pyplot.scatter(datax, datay, s=100, c=None, marker='o', cmap=None, norm=None,
                              #                         vmin=None, vmax=None, alpha=None, linewidths=None,
                              #verts=None, edgecolors=None, hold=None, data=None,
#**kwargs)

    plt.savefig("foo.png")

    """
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        
        #…(省略)…
        
        mappable0 = ax1.pcolormesh(X,Y,z, cmap='coolwarm', norm=LogNorm(vmin=1e1, vmax=1e3)) # ここがポイント！
        pp = fig.colorbar(mappable0, ax=ax1, orientation="vertical")
        pp.set_clim(1e1,1e3)
        pp.set_label(“color bar“, fontname="Arial", fontsize=10)
    """
