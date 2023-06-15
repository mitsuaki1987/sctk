.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

Introduction
============

This manual explains about the first-principles program package
"Superconducting-Toolkit" based on the density functional theory for
superconductors. Superconducting-Toolkit constructs the gap equation

.. math::
   
   \begin{align}
   \Delta_{n {\bf k}} = -\frac{1}{2} \sum_{n' {\bf k}'}
   \frac{K^{el}_{n {\bf k} n' {\bf k}'} + K^{\rm el-ph}_{n {\bf k} n' {\bf k}'}}{Z_{n {\bf k}}}
   \frac{\Delta_{n' {\bf k}'}}{E_{n' {\bf k}'}}
   \tanh\left( \frac{\beta E_{n' {\bf k}'}}{2} \right)
   \end{align}

from QuantumESPRESSO_ outputs and obtains the superconducting gap
function :math:`\Delta_{n {\bf k}}` by solving this equation. This program is
developed by Mitsuaki Kawamura in ISSP, Japan.
