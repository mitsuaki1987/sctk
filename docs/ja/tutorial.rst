.. _FermiSurfer: http://fermisurfer.osdn.jp/
.. _pw.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PW.html
.. _ph.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PH.html
.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

チュートリアル
==============

以下のチュートリアルは ``SCTK/examples/Al/`` 内で行う.

.. _scf:

電荷密度のSCF計算
-----------------

**入力ファイル**: scf.in

**プログラム**: pw.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ export OMP_NUM_THREADS=1
   $ mpiexec -np 29 PATH/pw.x -nk 29 -in scf.in > scf.out
        
**重要なパラメーター**

    calculation = "scf"
        pw.x_ でself-consistent計算を実行する.

.. note::

   次の :ref:`phonon`, :ref:`coulomb` はどちらを先にやっても変わらない.

.. _phonon:

フォノンおよび電子-フォノン相互作用の計算
-----------------------------------------

.. _ph:

フォノン振動数, 変移ポテンシャルの計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**入力ファイル**: ph.in

**プログラム**: ph.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ mpiexec -np 29 PATH/ph.x -nk 29 -in ph.in > ph.out

**重要なパラメーター**

    fildvscf = 'dv'
        変移ポテンシャルのファイル名. 次のステップの計算で使われる.

    ldisp = .true.
        一様 :math:`{\bf q}` グリッド上でフォノンを計算する.

    lshift_q = .true.
        :math:`\Gamma` 点の特異性を避けるために, :math:`{\bf q}` グリッドをずらす.

    nq1, nq2, nq3
        :math:`{\bf q}` グリッドの分割数. 

.. _elph:

電子-フォノン相互作用の計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**入力ファイル**: epmat.in

**プログラム**: ph.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ mpiexec -np 8 PATH/ph.x -nk 8 -in epmat.in > epmat.out

**重要なパラメーター**

    electron_phonon = "scdft_input"
        すでに計算されている変移ポテンシャル,
        ダイナミカルマトリックスを読み込み, 電子-フォノンバーテックスを計算する.

    elph_nbnd_min, elph_nbnd_max
        ギャップ方程式にはFermi準位近傍の電子-フォノン相互作用しか寄与しないので,
        計算コストを下げるために電子-フォノン相互作用を求めるバンドの本数を絞る.
        Fermi準位を含むバンドは, QuantumESPRESSO_ にある,
        `fermi_velocity.x <https://www.quantum-espresso.org/Doc/pp_user_guide/>`_
        によって調べることが出来る.

.. _coulomb:
   
遮蔽Coulomb/スピン揺らぎ相互作用の計算
--------------------------------------

.. _dense:

細かいkグリッドでのnon SCF計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**入力ファイル**: nscf.in

**プログラム**: pw.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ mpiexec -np 32 PATH/pw.x -nk 32 -in nscf.in > nscf.out
        
**重要なパラメーター**

   calculation = "nscf"
       Non Self-consistent 計算を行う.
   
   la2f = .true.
       Kohn-Shamエネルギーの情報を含むファイル pwscf.a2Fsave を出力する. 

   nbnd
       分極関数を計算するため, 非占有バンドも計算しておく必要がある.
       ただし, 半導体の計算の時ほどたくさんとる必要はない.
       目安は占有バンドと同程度である.

.. _twin:

遮蔽Coulomb相互作用計算のための波動関数の計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**入力ファイル**: twin.in

**プログラム**: pw.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ bash PATH/twingrid.x 4 4 4 >> twin.in
   $ mpiexec -np 32 PATH/pw.x -nk 32 -in twin.in > twin.out
        
**重要なパラメーター**

   calculation = "bands"
       この時 :math:`{\bf k}` 点メッシュに関しては上記のように
       :ref:`twingrid` の出力をファイル末尾にリダイレクトする.
       この時の, :math:`{\bf k}` 点メッシュは :ref:`elph` の ph.x_ のインプットの nq1,
       nq2, nq3 と同じにする. このインプットで pw.x_ を実行する.

遮蔽Coulomb相互作用の計算
~~~~~~~~~~~~~~~~~~~~~~~~~

**入力ファイル**: sctk.in

**プログラム**: :ref:`sctk.x <sctk>`

.. code-block:: bash

   $ mpiexec -np 32 PATH/sctk.x -nk 32 -in sctk.in > kel.out

**重要なパラメーター**

    :ref:`calculation = "kel" <kel>`
         遮蔽Coulom/スピン揺らぎ媒介相互作用の計算を行う.

    nq1, nq2,  nq3
         これらはこの前のステップの :math:`{\bf k}` グリッドと同じにする.
    
.. _scdftscf:
   
転移温度計算
------------

:ref:`phonon`, :ref:`coulomb` の計算が全て終了した段階でSCDFT計算を行う.

**入力ファイル**: sctk.in (一部書き換える)

**プログラム**: :ref:`sctk.x <sctk>`

.. code-block:: bash

   $ export OMP_NUM_THREADS=32
   $ mpiexec -np 1 PATH/sctk.x < sctk.in > tc.out
        
**重要なパラメーター**

    :ref:`calculation = "scdft_tc" <scdfttc>`
         この部分を書き換えて計算の種類を変える.
         ここでは2分法で :math:`T_c` を求める.

その他の解析
------------

パラメーター calculation を変えるといくつかの解析が行える.

**入力ファイル**: sctk.in (一部書き換える)

**プログラム**: :ref:`sctk.x <sctk>`

.. code-block:: bash

   $ mpiexec -np 1 PATH/sctk.x < sctk.in
        
**重要なパラメーター**

    :ref:`calculation = "lambda_mu_k" <lambdamuk>`
        FermiSurfer_
        でプロット可能な :math:`\lambda_{n {\bf {\bf k}}}` データを出力する.

    :ref:`calculation = "scdft" <scdft>`
        ある温度でのSCDFT計算を行い, 以下のポスト処理で必要なファイル delta.dat を出力する.

    temp
       温度(単位ケルビン). これを 0 もしくは負の値にすると,
       ゼロケルビン用の特別なアルゴリズムによる計算を行う.
       
    :ref:`calculation = "deltaf" <deltaf>`
        FermiSurfer_
        でプロット可能な :math:`\Delta_{n {\bf k}}` データを出力する.

    :ref:`calculation = "qpdos" <qpdos>`
         超伝導準粒子DOS. 計算時間は比較的長い.

.. note::

   SCTK/examples/MgB2 には別のチュートリアルがある.
   ただし, これらのチュートリアルは :math:`{\bf k}` 点数やバンド数, 擬ポテンシャルなどの精度が
   十分ではないことに注意.
