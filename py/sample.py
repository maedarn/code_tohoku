"""
参考:https://kyotogeopython.zawawahoge.com
"""

#出力
print("Hello, World!")

a = 3
b = 2*a
print(b)
print(a+b)

#型
a = 1 #int
b = 1.0 #float
c = "Spam!"
print(type(a))
print(type(b))
print(type(c))

#演算

x = 3
y = 2

print(x + y) # 足し算
print(x - y) # 引き算
print(x * y) # 掛け算
print(x / y) # 割り算 （Python 3から、intわるintでもfloatで返してくれるようになった。Python 2ではFortranなどと同じ）
print(x // y) # 割り算（整数部）
print(x % y) # 余り
print(x ** y) # べき乗

#理論式
x = 3
y = 2

print(x == y) # 等しい  x,yが数以外のときも用いることができる
print(x > y) # 大なり
print(x < y) # 小なり
print(x >= y) # 大なりイコール
print(x <= y) # 小なりイコール
print(x != y) # ノットイコール

a = (1 < 2) # a = True
b = (3**2 < 2**3) # b = False
print(a and b)
print(a or b)
print(a and (not b))

#配列
 #list []  -> my_list = [1.0, "orange", [1,2,3]] 自由
my_list = [1,2,3,4,5]
my_str = "Hello, world!"
print(len(my_list)) #len=ながさ
print(len(my_str))

print(my_list)
my_list.append(3.14) #最後に加える
print(my_list)
my_list.pop(1) # インデックスは0から始まることに注意,1の要素を消す
print(my_list)

 #tuple () , 書き換え不可
my_tuple = (1.0, [1,2,3], 3.14)
 #str（ストリング、文字列型）, " " or ' '
my_str = "abc"
new_str = my_str + "def" #後ろにつけられる

 #部分配列
my_list = [0,1,2,3,4,5,6,7,8,9]
slice1 = my_list[2:5]    # my_list[5]を含まないことに注目, 2->4
slice2 = my_list[1:8:3]  # 1から8まで（8を含まない）、増分3（2つおき）
slice3 = my_list[8:]     # 8から終わりまで
slice4 = my_list[::-1]   # 9から0まで1つずつ減っていく: stepが負の時はstartおよびstopのデフォルト値が入れ替わる

 #dictionary(整数ではなく文字列でインデクシングを行う),set(ダブった要素を同一視)
my_dict = {"Jim": 10, "Alex": 20, "Mike": 30} # dictionary
print(my_dict["Alex"]) #->20
my_set = {4,7,8,9,4,5,8,9,0} # set
print(my_set) #->{0, 4, 5, 7, 8, 9}

 #ユーザー定義型 (ex) NumPy ndarray : 演算なども自由に定義
import numpy as np       # NumPyモジュールのインポート: モジュールのインポートについては後述
my_arr1 = np.arange(10)   # NumPy Arrayの生成
my_arr2 = np.arange(0,20,2)
print(type(my_arr1))
print(my_arr1)
print(my_arr2)
print(my_arr1 + my_arr2)  # 要素ごとの足し算
print(my_arr1.sum())      # 和をとるメソッド

#型変換
print(str(123))            # int -> str
print(complex(-10.0))      # float -> complex
print(tuple(['a',1, 1.5])) # list -> tuple

#条件分岐: if, elif, else
"""
f [論理型]:
    [処理]         # ifの条件がTrueの場合
elif [論理型]: 
    [処理]         # ifの条件がFalseで、elifの条件がTrueの場合
else:
    [処理]         # ifもelifもFalseの場合
"""

time_now = 15

if (time_now < 7) or (23 < time_now): # False  f [論理型]:
    print("Sleeping")                 #  インデント(tab)    [処理]          ifの条件がTrueの場合,(スペース4つがインデント)
elif time_now < 12:                   # False  elif [論理型]:
    print("Good morning!")
elif time_now < 17:                   # True   else:
    print("Good afternoon!")
else:
    print("Good evening!")

#繰り返し: for, while
"""
for n in [配列型]:
    [処理]
while [論理型]:
    [処理]

range(start,end,interval)   # startからendまで（endは含まない）、間隔はinterval
range(start,end)            # 省略するとinterval = 1
range(end)                  # start = 0, interval = 1
"""
for i in range(10):
    cube = i**3
    print(cube)

#enumerate?

i = 10.0
while i > 1:
    print(i)
    i = i / 2

#関数:Pythonでは、Fortranにおけるsubroutineとfunctionの差はなく、まとめてfunction（関数）として扱われます。関数は基本的に値渡しです。
"""
def [関数名](引数):
    [処理]
    return [戻り値]
"""
 # 引数を2乗して返す関数
def square(x):
    return x**2

y = square(3) + square(4) # 3**2 + 4**2 = 25
print(y)

def bark():
    print("Bow-wow!")

bark()


 #注意 listの引き渡し
  # lis には参照が入る
lis = [10]

  # 代入文でも参照が渡される
newlis = lis

  # newlis を変更すると、その内容は lis にも伝播する
newlis[0] = 11
print("newlis", newlis)
print("lis", lis)

  # 配列の中身を変更する関数
def plusone(x):
    x[0] = x[0] + 1

  # この関数を呼ぶと、関数の中で行われた変更は、呼び出し元の変数 lis に伝播する
plusone(lis)
print("plusone", lis)

 #関数内で変数を書き換える場合、その影響は関数内にとどまります（ローカルスコープ）。関数外で定義された変数を書き換えたい場合は、global文を用いてその変数がグローバル変数であることを宣言する必要があります。
x = 10

def try_to_rewrite_1():
    global x
    x = x + 10

def try_to_rewrite_2():
    x = x + 10

print(x)
try_to_rewrite_1() # successful
print(x)
try_to_rewrite_2() # raises Error

#モジュール : 汎用的な関数や定数群などを異なる複数のスクリプトで用いるときは、モジュールを用いるのが便利です。 モジュールの作り方は簡単で、基本的には独立したスクリプト(my_module.pyとします)の中に関数や定数を順に書いていくだけです。

"""
module : my_module.pyを作った時の参照の仕方
モジュールで定義された変数や関数は、このスクリプト内の名前空間ではmy_module.[オブジェクト名]として扱われます。

pi = 3.141592

def area(r):
    return pi * r * r
"""
import my_module

radius = 5.0
S = my_module.area(radius) # my_module内のarea関数

pi = 0
print(pi)                  # このスクリプト内の定数pi
print(my_module.pi)        # my_module内の定数pi

import my_module as mod   # 別名をつける
print(mod.pi)

from my_module import pi  # モジュール内の特定のオブジェクトのみをインポート
print(pi)                  # この場合、モジュール内の選ばれたオブジェクトがこのスクリプトの名前空間内に定義される

 #同じdirectryにない時:importが参照するパスにmy_module.pyの場所を追加することができます(Pythonスクリプト上からsysモジュールを用いて次のようにパスを通すことができます)。
from sys import path  # pathは通っているパスが格納されているリストです
path.append("/home/username/where/module/is/located") # リストのappendメソッドで要素を追加

 #一度読まれたモジュールはimport文があっても再び読まれません。モジュールを修正して再び読み込んでほしい時には、セッションを再スタートする（IPythonなどを一度閉じる）か、importlibモジュール内のreload関数を用います。

import my_module
  # この時点でmy_module.pyを書き換えたとする
import my_module
  # ここではモジュールの更新が適用されていない
from importlib import reload
my_module = reload(my_module)
  # モジュールの更新が適用される

#コメント
# or """ """

#継続 -> \ ただしカッコ内では省ける


#matplotlib:matplotlibとは、Pythonで主に2次元の図を描画するときなどに使われる標準的なライブラリです。
"""
2通りの図の作成方法
実は、matplotlibで図を作成するのには、以下の2つの方法があります。

1. matplotlib.pyplotモジュールの関数を用いて図を作成する方法
plt.plot, plt.xlabelなど、pltモジュールの関数を用いて図の詳細を順々に指定してやる方法です。
非常に直感的でわかりやすい方法ですので、ちょっとデータを見たいときなどには重宝します。 
一方で、例えば「1つのファイルに複数の図を載せる」場合、plt.subplotで次の図の詳細を記述し始めると前の図の詳細の記述に戻れないなど、細かい図を作成するときにはやや不便です。

ここまで作成してきた図は、すべてこの1.の方法を用いて作成してきました。

2.matplotlib内のクラスを用いて図を作成する方法
matplotlib.pyplot.figure関数でFigureオブジェクトを作成し、図の詳細を追加していく方法です。 大まかな流れとしては、

①. matplotlib.pyplot.figure関数でFigureオブジェクトを生成する
②. Figureオブジェクトのメソッド(.add_subplot, .add_axesなど)を用いて、Axesオブジェクトを生成する
③. Axesオブジェクトのメソッド(.plot, .histなど)にデータを渡してプロットする
④. 必要ならば、Axesオブジェクトのメソッド(.set_xlim, set_titleなど)を用いて図を調整する

すなわち、1.の方法ではFigureオブジェクトが暗黙的に作られていたのに対し、2.の方法ではFigureオブジェクトを明示的に生成することになります。
記述量が多くなりますが、より自由に描写を操作することが可能になります。
"""


 # matplotlibライブラリのインポート
import matplotlib.pyplot as plt

 #↓この1行により、matplotlibで作成した図をJupyterに埋め込むことが出来ます
 #%matplotlib inline

 #一次元図
 import numpy as np
t=np.array([1,2,3,4])
y=t

plt.plot(t,y)
plt.show() #図を出力する

  #複数
t=np.array([1,2,3,4])
y1=t
y2=t**2
y3=t**3
plt.plot(t,y1)
plt.plot(t,y2)
plt.plot(t,y3)
plt.show()

  #その他ラベル替えなど
plt.figure(figsize=(6,8))      
plt.plot(t,y1,'o-',color='red',label='Primary',linewidth=5.0)
plt.plot(t,y2,'^--',color='green',label='Quadratic')
plt.plot(t,y3,'s',color='blue',label='Cubic')
plt.title('Various functions') 
plt.xlabel('t')                
plt.ylabel('y',fontsize=36)    
plt.xticks([1,3,4],[1,3,4.0],fontsize=24,rotation=45)             
plt.text(0.5,30.0,"rakugaki", fontsize=15, color="red") 
plt.legend(loc='upper left',fontsize=15)                
plt.xlim(0.0,5.0)              
plt.show()

 # ヒストグラム
x1 = np.random.randn(1000)  #平均0、標準偏差1の乱数を1000個生成する
x2 = np.random.randn(1000) + 2.0  #平均2、標準偏差1の乱数を1000個生成する

plt.hist(x1, normed=True, label = "x", bins = 20, range = (-3, 2), alpha = 0.5, color = "blue")
plt.hist(x2, normed=True, label = "y", bins = 20, range = (-3, 5), alpha = 0.5, color = "red")
plt.legend()
plt.show()
plt.show()


 # 方法２
 t=np.array([1,2,3,4])
y1=t
y2=t**2
y3=t**3

# Figureオブジェクトを生成
fig=plt.figure(figsize=(6,8))

# Axesオブジェクトを生成
ax=fig.add_subplot(1,1,1)

# 図のプロット
ax.plot(t,y1,'o-',color='red',label='Primary',linewidth=5.0)
ax.plot(t,y2,'^--',color='green',label='Quadratic')
ax.plot(t,y3,'s',color='blue',label='Cubic')

# 図の調整
ax.set_title('Various functions') 
ax.set_xlabel('t')                
ax.set_ylabel('y',fontsize=36)    
ax.set_xticks([1,3,4])
ax.set_xticklabels(['1','3','4.0'],fontsize=24,rotation=45)             
ax.text(0.5,30.0,"rakugaki", fontsize=15, color="red") 
ax.legend(loc='upper left',fontsize=15)                
ax.set_xlim(0.0,5.0)              

plt.show()



#Numpy 標準の演算が遅いというPythonの欠点を補うために開発された高速なベクトル演算をサポートするライブラリです。
"""
特徴
ベクトル演算をベースにした省メモリで高速な多次元配列（ndarrayオブジェクト）:Numpyで配列と言うとき、たいていndarray（N-dimensional array）オブジェクトのことを指します。
行列演算（逆行列など）や乱数生成、フーリエ変換を簡単に呼び出せる
PythonとC, C++, Fortranとの連携を可能にする低級インターフェイスの提供
バイナリやテキストでのファイル入出力を簡単にできる

注意
Numpyも万能ではなく、得意ではない処理もあります。それは、配列要素ごとに回すforループです。

Anacondaの場合「conda install numpy」
"""

 #インポート
import numpy as np    #慣習でnpと名付けることが多い
import matplotlib.pyplot as plt # 可視化用ライブラリ
 #一次元配列
 x = np.array([1., 0.01, 0.7])
 x
  #type()はPythonの組み込み関数で、引数が何のオブジェクトかを教えてくれます（dtypeとは違うものなので注意）。
  type(x)
  #shapeは、まさに配列の形状で、多次元配列であれば、軸に沿った要素の個数をそれぞれの軸について並べたタプルです。
  x.shape #xの行列の形
  #dtypeは、配列の型を表すものです
  x.dtype
  #ndimは、配列の次元を表します。
  x.ndim

 #多次元配列
  # ２×３の２次元配列を定義
  y = np.array([[1., 0.1, 3.],
              [3., 20., 0.07]])
  # 一行でも書いてよいですが、見やすい方法を選ぶといいでしょう。
  y = np.array([[1., 0.1, 3.], [3., 20., 0.07]])
  y
  # 複素数を入れた場合
  x = np.array([1.j, 3 + 6j, 9])
  x

  # 行0, 列1の値
  x[0, 1]
  # タプルと同じようにスライスを使える
  x[1:]
  # 多次元配列の場合
  x = np.arange(20).reshape(4,5)
  x
  # 列3, 4を抜き出す
  x[:, 3:]

 #演算
 x = np.array([1., 0.01, 0.7])
 x + 1 # ==> [ 2.  ,  1.01,  1.7 ]
 2 * x # ==> [ 2.  ,  0.02,  1.4 ]
 x / 1e5 # ==> [  1.00000000e-05,   1.00000000e-07,   7.00000000e-06]

 x = np.array([8., 0.3, 90])
 y = np.array([4., 0.2, 150])
 x + y # ==> [  12. ,    0.5,  240. ]
 x - y # ==> [  4. ,   0.1, -60. ]
 x * y # ==> [  3.20000000e+01,   6.00000000e-02,   1.35000000e+04]
 x / y # ==> [ 2. ,  1.5,  0.6]

 # ３×３行列を作成
 A = np.arange(9).reshape(3,3)
 A

 # ３×３の単位行列の２倍
 B = np.eye(3) * 2
 B

 # 要素ごとの積になる（行列積ではない！）
 A * B
 # 行列積
 A.dot(B)
 # Python3.5以降
 A @ B

#函数たち
 #np.array, np.asarray:上でも述べた関数です。リストやタプルを引数に取り、配列を返す関数です。
 x = np.array([1, 2, 5, 10])
 #np.arange:range()と似ていますが、さらに高機能な関数で、ndarrayを返します。
 np.arange(10) # 0から9まで
 np.arange(0, 10, 0.5) # 0から10まで0.5刻み（10は含まない）
 #np.ones, np.ones_like:np.onesは、引数で与えられたshapeで、全ての要素が１であるようなndarrayを返します。
 np.ones((2, 5)) # 多次元配列ver.
 np.ones_like(x) # xのようなshape, dtypeで、要素が全て１の配列
 #np.zeros, np.zeros_like:全ての要素を0にする他は、np.ones, np.ones_likeと同じです。
 #np.empty, np.empty_like:上の２つと似ていますが、こちらは初期化をせずに配列の領域だけを確保する関数です。巨大な配列を自分で初期化したい場合や、値は何でもいいが配列だけ作りたいときに使います。
 #np.linspace:ある範囲を均等に分割する関数です。
 np.linspace(0, 1, 10) # 0から1まで,3つ目の引数は個数を表す。
  #endpoint=Falseの場合。stopの値を含まない。
  np.linspace(0, 1, 10, endpoint=False)
 #np.eye, np.identity:この２つは単位行列を生成する関数です。
 #np.reshape:配列の次元を自在に変えられる関数です。全体のsizeを変更しないようにしなければならないこと、reshapeによって返ってくるndarrayは、元々のndarrayと同じメモリを参照している点に注意。
 b = np.arange(10) # オリジナル
 c = b.reshape(2,5) # reshapeされたビュー
 c[0, 0] = 100 # ビューを書き換える
 b # オリジナルが変わる

#ユニバーサル関数（ufunc）:Numpyでは、三角関数(np.sinなど）や指数関数（np.exp）のような関数はユニバーサル関数オブジェクト（略すとufunc）になっています。
 x = np.linspace(0.001, 5 * np.pi)
 plt.plot(x, np.sin(x) / x)
 """
# 平方根
np.sqrt([0, 2, 3, 5, 9])
array([0.        , 1.41421356, 1.73205081, 2.23606798, 3.        ])
# 自然対数
np.log([1, np.e, np.e**2])
array([0., 1., 2.])
# 底が10
np.log10([0.0001, 1.e-100, 1, 10, 100])
array([  -4., -100.,    0.,    1.,    2.])
# ラジアンから度への変換
np.rad2deg(np.linspace(0, np.pi, 10))
array([  0.,  20.,  40.,  60.,  80., 100., 120., 140., 160., 180.])
# 度からラジアンへの変換
np.deg2rad(np.linspace(0, 180, 10)) / np.pi
array([0.        , 0.11111111, 0.22222222, 0.33333333, 0.44444444,
       0.55555556, 0.66666667, 0.77777778, 0.88888889, 1.        ])

A = np.arange(6).reshape(3, 2)
A
array([[0, 1],
       [2, 3],
       [4, 5]])
A.sum(axis=0) # 縦方向の和
array([6, 9])
A.sum(axis=1) # 横方向の和
array([1, 5, 9])
A.sum() # 全要素の和
15

axisはsum,max,min,var,mean,stdなど多くの函数と併用
 """


 #ndarrayのdtypeのキャストはastypeを使います。
 a = np.arange(0, 2, 0.2)
 a.astype(np.int32) # ４バイト整数にキャスト

 #配列参照:コピーであるかビュー（参照）であるかを混乱することがあるので、チェックしておいてください。
  # スライスした部分配列（行1, 2）の要素に全て-1を代入する
  arr[1:3] = 0
  # オリジナルの配列
  arr

 # 一様乱数
 arr = np.random.rand(3, 4)
 arr
 # ３行目を書き換える
 arr[2] = arr[2] - 1
 arr
 # ２列めを書き換える
 arr[:, 1] = arr[:, 1] + arr[:, 2]
 arr
 #オリジナルの配列が変わってほしくないときは、copy関数で同じ配列を複製できます。
 arr_copied = arr.copy()
 arr_copied
 #代入時の注意点
  # ↓-100をaの全要素に代入するつもりで書いたら
  # 配列ですらなくなっている！？
  a = -100
  a
  # 正しくは
  a = np.arange(10)
  a[:] = -100
  a
 # ファンシーインデックス参照の一例（列0, 1を抜き出す）
 a[:, [0, 1]]
  # 列0, 1を抜き出す（コピーになる）
  b = a[:, [0, 1]]
  b

#C言語ライブラリのラッパー:他言語のライブラリの関数をPythonで呼び出せるようにしたものをラッパー（wrapper）といいます。
 #fft

 
#ファイル入出力
 #Numpyでテキストファイル入力を取り扱う関数として、loadtxtや、より高機能なgenfromtxtがあります。(pandas)
  # Anacondaで入れた場合はすでに使用可能。
  # エラーが出た人は、pipでインストールしてください。
  import pandas as pd

  #CSVファイルから読み取る
  file = """
  I1,I2,F1,I3,S1
  1,1,1.9029e+2,3,アツアツ
  3,2,8.24608e-1,9,萎え萎え
  """
  df = pd.read_csv(StringIO(file), delimiter=',')
  df
  # 実際は以下のようにパスを指定する。
  # df = pd.read_csv('data.csv', delimiter=',')

  #固定長幅テキストファイルから読み取る
  s = """         VALUE1            VALUE2         ID
  +1.87429510e+00   -4.44966444e-02   00030109
  -1.83643507e+00   -1.87712943e+00   04005971
  +5.11761193e-01   -5.57850439e-01   00000750
  """
  df = pd.read_fwf(StringIO(s), widths=(16, 19, 11))
  df
   # ndarrayが取得できるように
   df[['VALUE1', 'VALUE2']].values

  #テキストでの保存
   # ２次元配列
   x = np.linspace(0, 1, 10).reshape(5,2)
   # xをテキスト形式で保存
   np.savetxt("my_np.txt", x)

   """
   Fortranバイナリファイル読み込み
   ここでは、まず
   リトルエンディアンのヘッダなし4バイト浮動小数点バイナリ(俗に言うGrADS形式と呼ばれるもの)
   の中身を読み込んでみることにします。
   """
   import numpy as np
   N=10  #1レコード番号あたりに格納されているデータの数。
   M=20  #レコードの総数。
   f=open('filename.out','r')
   dty=np.dtype([('data','<'+str(N)+'f')])
   chunk=np.fromfile(f,dtype=dty,count=M)
   data=np.array([chunk[j]['data'] for j in range(0,M)])


#amiation
 # モジュールのインポート
 import matplotlib.animation as anm #<-これがanimationを作成する際に使用するモジュール。
 import matplotlib.pyplot as plt
 import numpy as np

 fig = plt.figure(figsize = (10, 6))
 x = np.arange(0, 10, 0.1)

 def update(i, fig_title, A):
     plt.cla()                      # 現在描写されているグラフを消去
     
     y = A * np.sin(x - i)
    plt.plot(x, y, "r")
    plt.title(fig_title + 'i=' + str(i))


    ani = anm.FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), \
                            interval = 100, frames = 32)
    # 作成されたアニメーションをファイルとして保存したいときには、以下の行も実行する。
    # gifファイルとして保存したいときには、予めImageMagickというアプリケーションをインストールした上で、writer = 'imagemagick' と指定する(UNIX系OSの場合)。
    ani.save("Sample.gif", writer = 'imagemagick')
    # mp4ファイルとして保存したいときには、予めffmpeg(動画と音声を変換することのできるUNIX系OS生まれのフリーソフトウェア)をインストールしておく(UNIX系OSの場合)。
    # ani.save("Sample.mp4")
