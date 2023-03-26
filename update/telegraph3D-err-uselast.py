# coding: UTF-8
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import h5py
import csv
import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation
from matplotlib.animation import ArtistAnimation

nx=128
ny=128
nz=128
nmesh=128
nd=19
nloop=50+1
step=40.0#/2.0/128.0

#dir1  = '/glv0/maedarn/test-telegraph/cg1-T20-mash256-L100-3pwe-4thcd/PHIINI/'
#dir2  = '/glv0/maedarn/test-telegraph/cg1-T10-mash256-L100-3pwe-4thcd/PHIINI/'
#dir3  = '/glv0/maedarn/test-telegraph/cg1-T50-mash256-L100-3pwe-4thcd/PHIINI/'
#dir4  = '/glv0/maedarn/test-telegraph/cg1-T5-mash256-L100-3pwe-4thcd/PHIINI/'
#dir5  = '/glv0/maedarn/test-telegraph/cg1-T100-mash256-L100-3pwe-4thcd/PHIINI/'

#dir1  = '/glv0/maedarn/test-telegraph/cg1-T20-mash256-L100-3pwe-4thmscl-test-t03-opr2/PHIINI/'
#dir1  = '/glv0/maedarn/test-telegraph/cg1-T20-mash512-L100-3pwe-4thcd-opr-2th/PHIINI/'
#dir1  = '/glv0/maedarn/test-telegraph/tel-mesh128-mv-sph-cg1-nit10/'
dir1  = '/glv0/maedarn/test-telegraph/tel-mesh128-kp-1mv/'

folder1="Phiwv"  #+ str(np.int(sample_frequency))
folder2="Phiexa"  #+ str(np.int(sample_frequency))

#folder1="Phiwv"  #+ str(np.int(sample_frequency))
#folder2="Phiexa"  #+ str(np.int(sample_frequency))

rms1 = [0] * nloop
rms2 = [0] * nloop
rms3 = [0] * nloop
rms4 = [0] * nloop
rms5 = [0] * nloop

rmsf = [0] * 5

rmsmd = [0] * 5

rmssk = [0] * 5
rmssk[0]=8.533865329809486866e-04
rmssk[1]=8.185834158211946487e-04
rmssk[2]=7.856506854295730591e-04
rmssk[3]=7.960871444083750248e-04
rmssk[4]=8.614736725576221943e-04
#nms = [ndx-4] * 5


rms6 = [0] * nloop
rms7 = [0] * nloop
rms8 = [0] * nloop
rms9 = [0] * nloop

ratio = [0] * nloop

b = np.arange(nloop)
b=b*step


h5file1last = h5py.File(dir1+'PHI128ITRt5T01000/NAllHDF'+"%03.f"%(nloop-1)+'.h5',"r")
h5file2last = h5py.File(dir1+'PHI128ITRt5T00500/NAllHDF'+"%03.f"%(nloop-1)+'.h5',"r")
h5file3last = h5py.File(dir1+'PHI128ITRt5T00200/NAllHDF'+"%03.f"%(nloop-1)+'.h5',"r")
h5file4last = h5py.File(dir1+'PHI128ITRt5T00030/NAllHDF'+"%03.f"%(nloop-1)+'.h5',"r")
h5file5last = h5py.File(dir1+'PHI128ITRt5T00050/NAllHDF'+"%03.f"%(nloop-1)+'.h5',"r")

Phiwv1to1last  = h5file1last[folder1+"/Phiwv1"].value
mx1last=np.max(Phiwv1to1last)
mn1last=np.min(Phiwv1to1last)

Phiwv1to2last  = h5file2last[folder1+"/Phiwv1"].value
mx2last=np.max(Phiwv1to2last)
mn2last=np.min(Phiwv1to2last)

Phiwv1to3last  = h5file3last[folder1+"/Phiwv1"].value
mx3last=np.max(Phiwv1to3last)
mn3last=np.min(Phiwv1to3last)

Phiwv1to4last  = h5file4last[folder1+"/Phiwv1"].value
mx4last=np.max(Phiwv1to4last)
mn4last=np.min(Phiwv1to4last)

Phiwv1to5last  = h5file5last[folder1+"/Phiwv1"].value
mx5last=np.max(Phiwv1to5last)
mn5last=np.min(Phiwv1to5last)



for i in range(1,nloop):
    
    h5file1 = h5py.File(dir1+'PHI128ITRt5T01000/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file2 = h5py.File(dir1+'PHI128ITRt5T00500/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file3 = h5py.File(dir1+'PHI128ITRt5T00200/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file4 = h5py.File(dir1+'PHI128ITRt5T00030/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file5 = h5py.File(dir1+'PHI128ITRt5T00050/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    '''
    h5file1 = h5py.File(dir1+'PHIT100/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file2 = h5py.File(dir1+'PHIT050/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file3 = h5py.File(dir1+'PHIT020/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file4 = h5py.File(dir1+'PHIT010/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file5 = h5py.File(dir1+'PHIT005/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    
    h5file6 = h5py.File(dir1+'PHIT015/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file7 = h5py.File(dir1+'PHIT030/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file8 = h5py.File(dir1+'PHIT106/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file9 = h5py.File(dir1+'PHIT107/NAllHDF'+"%03.f"%(i)+'.h5',"r")
    
    h5file6 = h5py.File(dir6+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file7 = h5py.File(dir7+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file8 = h5py.File(dir8+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file9 = h5py.File(dir9+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    '''
    
    #データ読み込み
    #Phiwv1to1  = h5file1[folder1+"/Phiwv1"].value
    Phiwv1to1  = h5file1[folder1+"/Phiwv1"].value
    Phiexa1 = h5file1[folder2+"/Phiexa"].value
    mx1=np.max(Phiexa1)
    mn1=np.min(Phiexa1)
    
    Phiwv1to2  = h5file2[folder1+"/Phiwv1"].value
    Phiexa2  = h5file2[folder2+"/Phiexa"].value
    mx2=np.max(Phiexa2)
    mn2=np.min(Phiexa2)
    
    Phiwv1to3  = h5file3[folder1+"/Phiwv1"].value
    Phiexa3  = h5file3[folder2+"/Phiexa"].value
    mx3=np.max(Phiexa3)
    mn3=np.min(Phiexa3)
    
    
    Phiwv1to4  = h5file4[folder1+"/Phiwv1"].value
    Phiexa4  = h5file4[folder2+"/Phiexa"].value
    mx4=np.max(Phiexa4)
    mn4=np.min(Phiexa4)
    
    Phiwv1to5  = h5file5[folder1+"/Phiwv1"].value
    Phiexa5  = h5file5[folder2+"/Phiexa"].value
    mx5=np.max(Phiexa5)
    mn5=np.min(Phiexa5)
    
    
    '''
    Phiwv1to6  = h5file6[folder1+"/Phiwv1"].value
    Phiexa6 = h5file6[folder2+"/Phiexa"].value
    mx6=np.max(Phiexa6)
    mn6=np.min(Phiexa6)
    
    Phiwv1to7  = h5file7[folder1+"/Phiwv1"].value
    Phiexa7  = h5file7[folder2+"/Phiexa"].value
    mx7=np.max(Phiexa7)
    mn7=np.min(Phiexa7)
    
    Phiwv1to8  = h5file8[folder1+"/Phiwv1"].value
    Phiexa8  = h5file8[folder2+"/Phiexa"].value
    mx8=np.max(Phiexa8)
    mn8=np.min(Phiexa8)
    
    
    Phiwv1to9  = h5file9[folder1+"/Phiwv1"].value
    Phiexa9  = h5file9[folder2+"/Phiexa"].value
    mx9=np.max(Phiexa9)
    mn9=np.min(Phiexa9)
    '''
    
    #spectrum = h5file[folder+"/spectrum"].value
    #data = np.loadtxt(dir53+'phi2D'+"%05.f"%(i)+'.dat', delimiter=',', unpack=True ,dtype='float')
    #rdata = np.reshape(data[a, :], (nd, ny ,nx))
    
    #err=0.25*(rdata[6,:,:]+rdata[7,:,:]+rdata[8,:,:]+rdata[9,:,:])-rdata[12,:,:]
    '''
    err1=Phiwv1to1- Phiexa1
    err2=Phiwv1to2- Phiexa2
    err3=Phiwv1to3- Phiexa3
    err4=Phiwv1to4- Phiexa4
    err5=Phiwv1to5- Phiexa5
    '''
    
    err1=Phiwv1to1- Phiwv1to1last
    err2=Phiwv1to2- Phiwv1to2last
    err3=Phiwv1to3- Phiwv1to3last
    err4=Phiwv1to4- Phiwv1to4last
    err5=Phiwv1to5- Phiwv1to5last
    
    '''
    err6=Phiwv1to6- Phiexa6
    err7=Phiwv1to7- Phiexa7
    err8=Phiwv1to8- Phiexa8
    err9=Phiwv1to9- Phiexa9
    '''
    rms1[i] = np.sqrt(np.mean(err1*err1))/(mx1-mn1)# -rmssk[0]
    rms2[i] = np.sqrt(np.mean(err2*err2))/(mx2-mn2)# -rmssk[1]
    rms3[i] = np.sqrt(np.mean(err3*err3))/(mx3-mn3)# -rmssk[2]
    rms4[i] = np.sqrt(np.mean(err4*err4))/(mx4-mn4)# -rmssk[3]
    rms5[i] = np.sqrt(np.mean(err5*err5))/(mx5-mn5)# -rmssk[4]
    '''
    rms1[i] = np.sqrt(np.mean(err1*err1))/(mx1last-mn1last)# -rmssk[0]
    rms2[i] = np.sqrt(np.mean(err2*err2))/(mx2last-mn2last)# -rmssk[1]
    rms3[i] = np.sqrt(np.mean(err3*err3))/(mx3last-mn3last)# -rmssk[2]
    rms4[i] = np.sqrt(np.mean(err4*err4))/(mx4last-mn4last)# -rmssk[3]
    rms5[i] = np.sqrt(np.mean(err5*err5))/(mx5last-mn5last)# -rmssk[4]
    '''
    '''
    rms6[i] = np.sqrt(np.mean(err6*err6))/(mx6-mn6)
    rms7[i] = np.sqrt(np.mean(err7*err7))/(mx7-mn7)
    rms8[i] = np.sqrt(np.mean(err8*err8))/(mx8-mn8)
    rms9[i] = np.sqrt(np.mean(err9*err9))/(mx9-mn9)
    '''
    #ratio[i]=rms1[i]/rms5[i]

'''
rms1[0]=100
rms2[0]=100
rms3[0]=100
rms4[0]=100
rms5[0]=100
'''
'''
rms1[0]=rms1[1]
rms2[0]=rms2[1]
rms3[0]=rms3[1]
rms4[0]=rms4[1]
rms5[0]=rms5[1]
'''
'''
rms6[0]=100
rms7[0]=100
rms8[0]=100
rms9[0]=100
'''
#rmsf[0] = rms1[99]
rmsf[0] = rms1[nloop-1]
rmsf[1] = rms2[nloop-1]
rmsf[2] = rms3[nloop-1]
rmsf[3] = rms4[nloop-1]
rmsf[4] = rms5[nloop-1]

#rmsmd[0] = rms1[0]#rms1[51]-rms1[50]
#rmsmd[1] = rms2[0]#rms2[51]-rms2[50]
#rmsmd[2] = rms3[0]#rms3[51]-rms3[50]
#rmsmd[3] = rms4[0]#rms4[51]-rms4[50]
#rmsmd[4] = rms5[0]#rms5[51]-rms5[50]

rmsmd[0] = 1.0#rms1[51]-rms1[50]
rmsmd[1] = 1.0#rms2[51]-rms2[50]
rmsmd[2] = 1.0#rms3[51]-rms3[50]
rmsmd[3] = 1.0#rms4[51]-rms4[50]
rmsmd[4] = 1.0#rms5[51]-rms5[50]
'''
rms6[0]=100i
rms7[0]=100
rms8[0]=100
rms9[0]=100

rms1[0]=rms1[1]
rms2[0]=rms2[1]
rms3[0]=rms3[1]
'''
ratio[0] = 1

fig = plt.figure()
#fig.subplots_adjust(bottom=0.21)
#fig.subplots_adjust(left=0.1)
#fig.subplots_adjust(right=0)
ax = fig.add_subplot(111)

plt.rcParams["font.family"] = "Times New Roman"      #全体のフォントを設定
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams['mathtext.default'] = 'it'
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams["xtick.major.size"] = 5                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 5                #y軸主目盛り線の長さ
plt.rcParams["font.size"] = 20                       #フォントの大きさ
plt.xlim([0.0, nloop*step])
#plt.xlim([1000.0, 1500]) # x方向の描画範囲を指定
plt.ylim([0.00001, 0.1]) # y方向の描画範囲を指定
#plt.rcParaims["legend.markerscale"] = 0.2
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.edgecolor"] = 'black'

#arr += -124
#rms1 += -rmssk[0]
#rms2 += -rmssk[1]
#rms3 += -rmssk[2]
#rms4 += -rmssk[3]
#rms5 += -rmssk[4]

plt.yscale('log')
#plt.plot(b,rms9, color='orange' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5 × 10^{-6}$")
#plt.plot(b,rms8, color='darkviolet' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5 × 10^{-5}$")

plt.plot(b,rms1, color='red' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=0.5$")
plt.plot(b,rms2, color='blue' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=1.0$")
#plt.plot(b,rms7, color='olive' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=1.667")
plt.plot(b,rms3, color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=2.5$") #2.5
#plt.plot(b,(rms1-rms1[50])/rmsmd[0], color='red' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=0.5$") #2.5
#plt.plot(b,rms6, color='grey' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=3.333")
plt.plot(b,rms4, color='black' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5.0$")
#plt.plot(b,rms3, color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="T50")
plt.plot(b,rms5, color='magenta' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=10$")
#plt.plot(b,rms1, color='red' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=0.5$")
#plt.plot(b,rms5, color='magenta' , linestyle = "solid", markersize=2.5, linewidth = 1,label="T100")
#plt.plot(b,rms1, color='red' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$v=0.02c_g$")
#plt.plot(b,(rms2-rms2[50])/rmsmd[1], color='blue' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=1.0$")
#plt.plot(b,(rms3-rms3[50])/rmsmd[2], color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=2.5$")
#plt.plot(b,(rms4-rms4[50])/rmsmd[3], color='black' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5.0$")
#plt.plot(b,(rms5-rms5[50])/rmsmd[4], color='magenta' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=10.0$")
#plt.plot(b,rms2, color='blue' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=1.0$")
#plt.plot(b,rms3, color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=2.5$")
#plt.plot(b,rms4, color='black' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5.0$")
#plt.plot(b,rms5, color='magenta' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=10.0$")

#plt.plot(b,ratio, color='blue' , linestyle = "solid", markersize=2.5, linewidth = 1,label="T10")


ax.set_xlabel('step')
ax.set_ylabel('err')

plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=9,frameon=False)

fig.tight_layout()
#グラフ表示
plt.show()
#plt.figure(figsize=(5, 6))

# 保存するPDFファイル名
#pp = PdfPages(dir1+'test1.pdf')
#pp = PdfPages('/home/maedarn/newkaiseki/test'+"%03.f"%(nmesh)+'-T20-r10-opr-2th2.pdf')
#pp = PdfPages('/home/maedarn/newkaiseki/tel-err064.pdf')
pp = PdfPages('/home/maedarn/newkaiseki/an-tel/tel-err-mv-itr-kp-cnt-last'+"%03.f"%(nmesh)+'.pdf')
# 画像をPDFとして保存する
pp.savefig(fig)
# PDFの保存終了
pp.close()
#plt.savefig("R-Vesc.png")

#ratio[100]

#rmsf[0]=np.mean(rms1)
#rmsf[1]=np.mean(rms2)
#rmsf[2]=np.mean(rms3)
#rmsf[3]=np.mean(rms4)
#rmsf[4]=np.mean(rms5)


np.savetxt('/home/maedarn/newkaiseki/an-tel/tel-err-mv-itr-kp-cnt-last'+"%03.f"%(nmesh)+'.txt', rms1)
#rms1[40],rms1[40]
