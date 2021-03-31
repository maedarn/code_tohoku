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

nx=64+4
ny=64+4
nz=64+4
nd=19
nloop=100+1
step=20

dir1  = '/glv0/maedarn/test-telegraph/cg1-T20-mash64-L100-3pwe-4thcd/PHIINI/'
dir2  = '/glv0/maedarn/test-telegraph/cg1-T10-mash64-L100-3pwe-4thcd/PHIINI/'
dir3  = '/glv0/maedarn/test-telegraph/cg1-T50-mash64-L100-3pwe-4thcd/PHIINI/'
dir4  = '/glv0/maedarn/test-telegraph/cg1-T5-mash64-L100-3pwe-4thcd/PHIINI/'
dir5  = '/glv0/maedarn/test-telegraph/cg1-T100-mash64-L100-3pwe-4thcd/PHIINI/'



folder1="Phiwv"  #+ str(np.int(sample_frequency))
folder2="Phiexa"  #+ str(np.int(sample_frequency))


rms1 = [0] * nloop
rms2 = [0] * nloop
rms3 = [0] * nloop
rms4 = [0] * nloop
rms5 = [0] * nloop

rmsf = [0] * 5

#nms = [ndx-4] * 5

'''
rms6 = [0] * nloop
rms7 = [0] * nloop
rms8 = [0] * nloop
rms9 = [0] * nloop
'''
ratio = [0] * nloop

b = np.arange(nloop)
b=b*step
for i in range(1,nloop):
    
    h5file1 = h5py.File(dir1+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file2 = h5py.File(dir2+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file3 = h5py.File(dir3+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file4 = h5py.File(dir4+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file5 = h5py.File(dir5+'NAllHDF'+"%03.f"%(i)+'.h5',"r")

    '''
    h5file6 = h5py.File(dir6+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file7 = h5py.File(dir7+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file8 = h5py.File(dir8+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    h5file9 = h5py.File(dir9+'NAllHDF'+"%03.f"%(i)+'.h5',"r")
    '''
    
    #データ読み込み
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
    err1=Phiwv1to1- Phiexa1
    err2=Phiwv1to2- Phiexa2
    err3=Phiwv1to3- Phiexa3
    err4=Phiwv1to4- Phiexa4
    err5=Phiwv1to5- Phiexa5
    
    '''
    err6=Phiwv1to6- Phiexa6
    err7=Phiwv1to7- Phiexa7
    err8=Phiwv1to8- Phiexa8
    err9=Phiwv1to9- Phiexa9
    '''
    rms1[i] = np.sqrt(np.mean(err1*err1))/(mx1-mn1)
    rms2[i] = np.sqrt(np.mean(err2*err2))/(mx2-mn2)
    rms3[i] = np.sqrt(np.mean(err3*err3))/(mx3-mn3)
    rms4[i] = np.sqrt(np.mean(err4*err4))/(mx4-mn4)
    rms5[i] = np.sqrt(np.mean(err5*err5))/(mx5-mn5)
    '''
    rms6[i] = np.sqrt(np.mean(err6*err6))/(mx6-mn6)
    rms7[i] = np.sqrt(np.mean(err7*err7))/(mx7-mn7)
    rms8[i] = np.sqrt(np.mean(err8*err8))/(mx8-mn8)
    rms9[i] = np.sqrt(np.mean(err9*err9))/(mx9-mn9)
    '''
    ratio[i]=rms1[i]/rms5[i]


rms1[0]=100
rms2[0]=100
rms3[0]=100
rms4[0]=100
rms5[0]=100

rmsf[0] = rms1[100]
rmsf[1] = rms2[100]
rmsf[2] = rms3[100]
rmsf[3] = rms4[100]
rmsf[4] = rms5[100]

'''
rms6[0]=100
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
plt.xlim([0.0, nloop*step]) # x方向の描画範囲を指定
plt.ylim([0.001, 0.1]) # y方向の描画範囲を指定
plt.rcParams["legend.markerscale"] = 2
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.edgecolor"] = 'black'


plt.yscale('log')
#plt.plot(b,rms9, color='orange' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5 × 10^{-6}$")
#plt.plot(b,rms8, color='darkviolet' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=5 × 10^{-5}$")

plt.plot(b,rms5, color='magenta' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa=0.5$")
plt.plot(b,rms3, color='olive' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=1.0")
#plt.plot(b,rms7, color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=1.667")
plt.plot(b,rms1, color='red' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=2.5")
#plt.plot(b,rms6, color='grey' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=3.333")
plt.plot(b,rms2, color='blue' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=5.0")
#plt.plot(b,rms3, color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="T50")
plt.plot(b,rms4, color='black' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$\kappa$=10")
#plt.plot(b,rms5, color='magenta' , linestyle = "solid", markersize=2.5, linewidth = 1,label="T100")
#plt.plot(b,rms1, color='red' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$v=0.02c_g$")
#plt.plot(b,rms2, color='blue' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$v=0.1c_g$")
#plt.plot(b,rms3, color='green' , linestyle = "solid", markersize=2.5, linewidth = 1,label="$v=0.5c_g$")

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
pp = PdfPages('/home/maedarn/newkaiseki/test64.pdf')
# 画像をPDFとして保存する
pp.savefig(fig)
# PDFの保存終了
pp.close()
#plt.savefig("R-Vesc.png")

ratio[100]

np.savetxt('/home/maedarn/newkaiseki/mesh64.txt', rmsf)
