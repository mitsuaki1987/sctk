.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

はじめに
========

この文書では超伝導密度汎関数理論(SCDFT)に基づく第一原理計算プログラムパッケージ
「Superconducting Toolkit」についての解説を行っている.
Superconducting Toolkitは QuantumESPRESSO_ のアウトプットを用いてギャップ方程式

.. math::
   
   \begin{align}
   \Delta_{n {\bf k}} = -\frac{1}{2} \sum_{n' {\bf k}'}
   \frac{K^{el}_{n {\bf k} n' {\bf k}'} + K^{\rm el-ph}_{n {\bf k} n' {\bf k}'}}{Z_{n {\bf k}}}
   \frac{\Delta_{n' {\bf k}'}}{E_{n' {\bf k}'}}
   \tanh\left( \frac{\beta E_{n' {\bf k}'}}{2} \right)
   \end{align}

を構成し, それを解いてギャップ関数 :math:`\Delta_{n {\bf k}}` を計算する.
このソフトウェアは東京大学物性研究所の河村光晶によって開発された.
