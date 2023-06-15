.. _pw.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PW.html
.. _ph.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PH.html

波数グリッドおよびバンド範囲の指定について
==========================================

本計算では様々な波数グリッドやバンドを指定する範囲が出てくるため
分かりにくい場合がある.
ここではそれらについて説明をする.

バンドの範囲
------------

- :ref:`scf`, :ref:`ph` でのバンド上限 ``nbnd``\ (scf)

     pw.x_ が価電子数から自動決定するものを使えばよい.
     したがって入力ファイルで指定する必要はない.

- :ref:`coulomb` でのバンド上限 ``nbnd``\ (:math:`K^{el}`)

     典型的には ``nbnd``\ (:math:`K^{el}`)を ``nbnd``\ (scf)の2倍程度にする.

     :ref:`sctk.x <sctk>` の計算コストは, ``nbnd``\ (:math:`K^{el}`)の2乗に比例する.

- :ref:`elph` での ``elph_nbnd_min, elph_nbnd_max``

     典型的にはFermi準位を含むバンドの下限と上限
     (`fermi_velocity.x <https://www.quantum-espresso.org/Doc/pp_user_guide/>`_
     によって調べることが出来る)だが,
     極めて大きいフォノンの振動数を持つ物質に対してはより広くとる必要がある.

- :ref:`scdfttc` における, 電子-電子Coulomb項のバンドの上限と下限 ``fbee, lbee``

     ``fbee=1, lbee=nbnd``\ (:math:`K^{el}`), すなわちデフォルトのままで良い.
     バンド数に対する収束を調べる場合などにそこから減らしてみる.

- :ref:`deltaf` で出力されるバンドの上限と下限 ``fbfs, lbfs``
   
     (バンド交差を無視した上での)Fermi面を横切るバンドの上限と下限が
     プログラムによって自動計算される.

これらのバンドの大小関係は次のようになる.

1 :math:`\leq` ``fbee`` :math:`\leq` ``elph_nbnd_min`` :math:`\leq`
``fbfs`` :math:`\leq` ``lbfs`` :math:`\leq` ``elph_nbnd_max``
:math:`\approx` ``nbnd``\ (scf) :math:`\leq` ``lbee`` :math:`\leq` ``nbnd``\ (:math:`K^{el}`) 

波数グリッド
------------

- :ref:`scf`, :ref:`ph` での, 電子状態に関する波数グリッド

   pw.x_ のインプットファイルで, 次のように指定される.
  
   ::
      
      K_POINTS automatic 
      {nk1} {nk2} {nk3} 0 0 0


   このグリッドで得られる :math:`{\bf k}` 点数 :math:`N_{\bf k}^{\rm smooth}` に対する各プログラムの
   計算コストの依存性は次のようになる.

   - :ref:`scf`, :ref:`ph` の計算コストは :math:`N_{\bf k}^{\rm smooth}` に比例する.
      
- :ref:`phonon`, :ref:`twin` での波数グリッド

   ph.x_ のインプット ``nq1, nq2, nq3`` と,
   :ref:`twingrid` の引数,
   および :ref:`elph` での ``nk1, nk2, nk3`` は同じにしなければならない.

   このグリッドで得られる :math:`{\bf q}` 点数 :math:`N_{\bf q}` に対する各プログラムの
   計算コストの依存性は次のようになる.

   - pw.x_ の計算コストは :math:`N_{\bf q}` に比例する.
      
   - :ref:`ph` での(すべての :math:`{\bf q}` での計算を合わせた)計算コストは :math:`N_{\bf q}` に比例する.

   - :ref:`elph` での(すべての :math:`{\bf q}` での計算を合わせた)計算コストは :math:`N_{\bf q}^2` に比例する.

   - :ref:`sctk.x <sctk>` での(すべての :math:`{\bf q}` での計算を合わせた)計算コストは :math:`N_{\bf q}^2` に比例する.

- :ref:`dense` での波数グリッド :ref:`[1] <ref>`

   これは状態密度を計算するときと同程度に細かい :math:`{\bf k}` グリッドを取る必要がある.
   このグリッドで得られる :math:`{\bf k}` 点数 :math:`N_{\bf k}^{\rm dense}` に対する各プログラムの
   計算コストの依存性は次のようになる.

   - :ref:`scdft` および :ref:`sctk.x <sctk>` の計算コストは :math:`N_{\bf k}^{\rm dense}`
     にあまり影響されない.

   - :ref:`deltaf` の計算コストは :math:`N_{\bf k}^{\rm dense}` に比例する.

これらの波数グリッドの大小関係は次のようになる.

:math:`N_{\bf q} \leq N_{\bf k}^{\rm smooth} \leq N_{\bf k}^{\rm dense}`
