.. _FermiSurfer: http://fermisurfer.osdn.jp/
.. _pw.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PW.html
.. _ph.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PH.html

.. _sctk:

プログラムの説明
================

Superconducting Toolkitのメインプログラム sctk.x の機能は次のとおりである。

-  前処理として,
   QEで計算したKohn-Sham軌道から 電子間Coulom相互作用行列
   :math:`K^{el}_{n {\bf k} n' {\bf k}'}` を計算する.
-  ギャップ方程式を解いて超伝導ギャップ :math:`\Delta_{n {\bf k}}` を計算する.
-  ギャップ方程式を解いて超伝導ギャップ :math:`\Delta_{n {\bf k}}` を計算し、２分法で転移温度を計算する.
-  後処理として、ギャップ方程式を非自己無撞着に計算し, FermiSurfer_
   でプロット可能なファイルを出力する.
-  後処理として、ギャップ方程式を非自己無撞着に計算してから, 準粒子DOSを計算する.
-  後処理として、ギャップ方程式を非自己無撞着に計算してから, 超音波吸収係数を計算する.

この他に、 :ref:`Coulomb相互作用の計算 <kel>` で用いるKohn-Sham軌道を計算する
2重 :math:`{\bf k}` 点グリッドを生成するためのスクリプト :ref:`twingrid` がある。

使い方
------

.. code-block:: bash

    mpiexec options PATH/sctk.x -in input_file
        
入力ファイル書式
----------------

入力ファイルの書式は次のとおりである.

::

    &CONTROL
      prefix = 
      outdir = 
      calculation = 
    /
    &KEL
      start_q =
      last_q =
      nci =
      laddxc =
      ecutfock =
      nq1 =
      nq2 =
      nq3 =
      lsf
    /
    &SCDFT
      temp =
      fbee =
      lbee =
      xic =
      nmf =
      nx =
      ne =
      emin =
      emax =
      electron_maxstep =
      conv_thr =
      fildyn =
      spin_fluc =
    /
        
これらのパラメーターはネームリスト内では順不同に指定できる.
また省略した場合はデフォルト値が使われる.

ネームリスト CONTROL
~~~~~~~~~~~~~~~~~~~~

============ ========= ============ ===================================================================
パラメーター 型        デフォルト値 説明
============ ========= ============ ===================================================================
prefix       文字列    "pwscf"      pw.x_ の入力の ``prefix`` と同じにする.
outdir       文字列    "./"         pw.x_ の入力の ``outdir`` と同じにする.
calculation  文字列    "kel"        計算の種類
============ ========= ============ ===================================================================

ネームリスト KEL
~~~~~~~~~~~~~~~~

``calculation = "kel"`` としたときに使われるパラメーター。

============= ========= =============== ===================================================================
パラメーター  型        デフォルト値    説明
============= ========= =============== ===================================================================
start_q       正の整数  1               遮蔽Coulomb/スピンゆらぎ相互作用を計算する :math:`{\bf q}` 点の始点
last_q        正の整数  既約            遮蔽Coulomb/スピンゆらぎ相互作用を計算する :math:`{\bf q}` 点の終点
                        :math:`{\bf q}`
                        点数
laddxc        0 or 1    0               遮蔽を計算するときの近似のレベルを指定する. 0 : RPA, 1: 断熱的LDA.
lsf           0 or 1    0               スピンゆらぎ相互作用 :ref:`[3] <ref>` を計算する(1)かしない(0)か
ecutfock      実数      pw.x_ での値    分極関数を計算するときの平面波カットオフ [Ry].
nq1, nq2, nq3 正の整数  a2Fsaveの       波動関数データの :math:`{\bf k}` 点メッシュ.
                        :math:`{\bf k}` :ref:`twingrid` の入力と同じにしなければならない. 
                        点メッシュ数と
                        同じ        
nci           正の整数  5               遮蔽Coulomb相互作用を計算する松原振動数の数.
============= ========= =============== ===================================================================

ネームリスト SCDFT
~~~~~~~~~~~~~~~~~~

``calculation = "scdft"`` としたときなどに使われるパラメーター。

================ ======== ============ ===================================================================
パラメーター     型       デフォルト値 説明
================ ======== ============ ===================================================================
temp             正の実数 0.1          温度. 単位ケルビン.
fbee             正の整数 1            全バンドのうち, ギャップ方程式の計算に含める
                                       一番初めのバンド.
lbee             正の整数 pw.x_ のnbnd ギャップ方程式の計算に含める最後のバンド.
xic              実数     -1.0         ギャップ関数外挿法に用いるパラメーター. 単位 Ry.
                                       これを ``0.0`` 未満にするとギャップ関数外挿法を使わない.
                                       デフォルトでは外挿法を使わない設定になっている.
nmf              整数     10           Comlombカーネルの計算で用いる松原振動数積分に用いる点の数.
                                       ``0`` にすると静的なCoulomb相互作用のみをつかう. また,
                                       負の値にするとCoulomb相互作用項を0として
                                       (フォノン項のみを考慮して)計算する.
nx               正の整数 100          フェルミ面近傍のバンドの付加的エネルギーグリッドのグリッド数.
ne               正の整数 50           :ref:`準粒子DOS計算 <qpdos>` のみで使用.
                                       準粒子DOSを計算するエネルギー点数.
emin             正の実数 1.0e-7       フェルミ面近傍のバンドの付加的エネルギーグリッドのためのパラメータ.
                                       単位 Ry.
emax             正の実数 5.0          :ref:`準粒子DOS計算 <qpdos>` のみで使用.
                                       準粒子DOSを計算するエネルギーグリッドの上限. 単位 meV.
electron_maxstep 正の整数 100          ギャップ方程式を反復法で解くときの反復回数の上限数.
conv_thr         正の実数 1.0e-15      ギャップ方程式を反復法で解くときの,
                                       新旧のギャップ関数の差の2乗平均に対する収束判定のしきい値. 単位 Ry.
filedyn          文字列   "matdyn"     ph.x_ の filedyn と同じにしなければならない。
spin_fluc        論理型   .False.      .True. にするとスピン揺らぎ :ref:`[3] <ref>` を含める。
scdft_kernel     正の整数 1            1: Lüders2005 :ref:`[4] <ref>`, 2: Sanna2020 :ref:`[5] <ref>`
lz_coulomb       論理型   .False.      Coulomb renormalization :ref:`[6] <ref>`
================ ======== ============ ===================================================================

入出力ファイル
--------------

sctk.xに関連するファイルは次の通りである。

.. _xml:

{prefix}.xml
~~~~~~~~~~~~

格子定数等の情報を含む. pw.x_ により生成される.

.. _a2fsave:

{prefix}.a2Fsave
~~~~~~~~~~~~~~~~

通常のDFT計算で求めたKohn-Shamエネルギーやその :math:`{\bf k}` メッシュ情報,
対称操作を含む. pw.x_ で la2f=.true. とすると生成される.

.. _wfc:

{prefix}.save/wfc\*.dat
~~~~~~~~~~~~~~~~~~~~~~~

各 :math:`{\bf k}` 点のKohn-Sham軌道.
\* には :math:`{\bf k}` 点の番号が入る.
pw.x_ により生成される.

.. _veldat:

vel\*.dat
~~~~~~~~~

各 :math:`{\bf q}` 点での遮蔽Coulomb相互作用のChebyshev補間の係数.
\* には :math:`{\bf q}` 点の番号が入る.
sctk.x で :ref:`calculation="kel" <kel>` とすると出力される.

.. _elphdat:

elph\*.dat
~~~~~~~~~~

各 :math:`{\bf q}` 点での電子-フォノン相互作用, フォノン振動数.
\* には :math:`{\bf q}` 点の番号が入る.
ph.x_ で electron_phonon="scdft_input" とすると作られる.

.. _lambdafrmsf:

lambda.frmsf, mu.frmsf
~~~~~~~~~~~~~~~~~~~~~~

くりこみ因子 :math:`\lambda_{n {\bf k}}` のFermi面上での値をプロットするための,
FermiSurfer_ 用データファイル.
sctk.x で :ref:`calculation="lambda_mu_k" <lambdamuk>` とすると出力される.

.. _deltadat:

delta.dat
~~~~~~~~~

超伝導ギャップ関数 :math:`\Delta_{n {\bf k}}`.
対応するKohn-Shamエネルギー :math:`\xi_{n {\bf k}}`, 積分重み, バンド番号,
:math:`{\bf k}` 点番号, 繰りこみ因子 :math:`Z_{n {\bf k}}`
sctk.x で :ref:`calculation="scdft" <scdft>` とすると出力される.

.. _qpdosdat:

qpdos.dat
~~~~~~~~~

第1列:準粒子エネルギー(単位 meV),
第2列:準粒子状態密度(単位 Ry\ :math:`^{-1}`).
sctk.x で :ref:`calculation="qpdos" <qpdos>` とすると出力される.

.. _deltafrmsf:

delta.frmsf, Z.frmsf
~~~~~~~~~~~~~~~~~~~~

超伝導ギャップ関数 :math:`\Delta_{n {\bf k}}` およびくりこみ因子 :math:`Z_{n {\bf k}}`
のFermi面上での値をプロットするための,
FermiSurfer_ 用データファイル.
sctk.x で :ref:`calculation="deltaf" <deltaf>` とすると出力される.

計算の種類
----------

パラメーター calculation に次の文字列を入れて、計算の種類を指定する。

.. _kel:

kel : 遮蔽Coulomb/スピン揺らぎ媒介相互作用
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

プログラムを実行しているディレクトリ内に,
次のものを用意しておく必要がある.

-  :ref:`xml`
-  :ref:`a2fsave`
-  :ref:`wfc`

プログラムを実行したディレクトリに, 次のものが作られる.

-  :ref:`veldat`

.. _lambdamuk:

lambda_mu_k : 軌道ごとの電子フォノンパラメーター
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

プログラムを実行しているディレクトリ内に,
次のものを用意しておく必要がある.

-  :ref:`xml`
-  :ref:`a2fsave`
-  :ref:`elphdat`   
-  :ref:`veldat`

プログラムを実行したディレクトリに, 次のものが作られる.

-  :ref:`lambdafrmsf`

.. _scdft:

scdft : ある温度でのSCDFT計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

プログラムを実行しているディレクトリ内に,
次のものを用意しておく必要がある.

-  :ref:`xml`
-  :ref:`a2fsave`
-  :ref:`elphdat`   
-  :ref:`veldat`

プログラムを実行したディレクトリに, 次のものが作られる.

-  :ref:`deltadat`

.. _scdfttc:

scdft_tc : ２分法による転移温度の自動計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

プログラムを実行しているディレクトリ内に,
次のものを用意しておく必要がある.

-  :ref:`xml`
-  :ref:`a2fsave`
-  :ref:`elphdat`   
-  :ref:`veldat`

プログラムを実行したディレクトリに, 次のものが作られる.

-  :ref:`deltadat`

.. _qpdos:

qpdos : 準粒子状態密度
~~~~~~~~~~~~~~~~~~~~~~

プログラムを実行しているディレクトリ内に,
次のものを用意しておく必要がある.

-  :ref:`xml`
-  :ref:`a2fsave`
-  :ref:`elphdat`   
-  :ref:`veldat`
-  :ref:`deltadat`

プログラムを実行したディレクトリに, 次のものが作られる.

-  :ref:`qpdosdat`

.. _deltaf:

deltaf : Fermi面上でのギャップ関数を計算しFermiSurfer用ファイルを出力する
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

プログラムを実行しているディレクトリ内に,
次のものを用意しておく必要がある.

-  :ref:`xml`
-  :ref:`a2fsave`
-  :ref:`elphdat`   
-  :ref:`veldat`
-  :ref:`deltadat`

プログラムを実行したディレクトリに, 次のものが作られる.

-  :ref:`deltafrmsf`

ultrasonic : 超音波吸収係数
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _twingrid:
   
twingrid.x
----------

sctk.x において :ref:`calculation="kel" <kel>` でCoulomb相互作用を計算するときのKohn-Sham軌道の計算(pw.x_)
において用いる2重 :math:`{\bf k}` グリッドを生成するスクリプト.

使い方
~~~~~~

.. code-block:: bash

    $ bash PATH/twingrid.x nk1 nk2 nk3 >> input_file_for_pw
        

*nk1*, *nk2*, *nk3* はそれぞれの逆格子ベクトルの方向の :math:`{\bf k}` 点分割数.

標準出力
~~~~~~~~

次のように標準出力される.

::

    K_POINTS crystal
    Total_number_of_k
    k_vector1 1.0
    k_vector2 1.0
    k_vector3 1.0
     :
        
これにより,  :math:`\Gamma` 点を含むグリッドと,
そこから半グリッドぶんずらしたグリッド上の :math:`{\bf k}` 点がセットで生成される.
上の使い方では,
この標準出力を pw.x_ の入力ファイルの末尾にリダイレクトしている.


