.. _pw.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PW.html
.. _ph.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PH.html

Setting of wave-number grid and band range
==========================================

In the calculation with this package,
there are many kind of the wave-number grid and the band range;
It may confuse us.
In this chapter, the relation and the difference between them
are described.

Range of bands
--------------

- The upper limit in :ref:`scf`, :ref:`ph` : ``nbnd``\ (scf)

     We should use the number specified automatically by ``pw.x``.
     Therefore, we do not have to write explicitly in the input file.

- The upper limit in :ref:`coulomb`: ``nbnd``\ (:math:`K^{el}`)

     Typically, ``nbnd``\ (:math:`K^{el}`) should be roughly the double of ``nbnd``\ (scf).

     The numerical cost of :ref:`sctk.x <sctk>`  is proportional to the square of ``nbnd``\ (:math:`K^{el}`).

- ``elph_nbnd_min, elph_nbnd_max`` in :ref:`elph`

     In almost cases, they are equal to the lower- and the upper limit of bands
     that contain the Fermi level 
     (These limit can be obtained by
     `fermi_velocity.x <https://www.quantum-espresso.org/Doc/pp_user_guide/>`_ ).
     For materials that have extremely large phonon frequencies,
     this band range must be wider than ordinary cases.

- The lower- and the upper- limit for the electron-electron Coulomb term in :ref:`scdfttc` : ``fbee, lbee``

     The default value [``fbee=1, lbee=nbnd``\ (:math:`K^{el}`)] is recommended.
     When we check the convergence about the number of :math:`{\bf k}` point,
     we reduce them from the default value.

- The lower- and the upper limit of bands printed by :ref:`deltaf` : ``fbfs, lbfs``
   
     They are the lower- and the upper limit of bands
     that contain the Fermi level (with non-crossing approximation).
     They are computed automatically.

The relation of magnitude of these bands becomes as follows:

1 :math:`\leq` ``fbee`` :math:`\leq` ``elph_nbnd_min`` :math:`\leq`
``fbfs`` :math:`\leq` ``lbfs`` :math:`\leq` ``elph_nbnd_max``
:math:`\approx` ``nbnd``\ (scf) :math:`\leq` ``lbee`` :math:`\leq` ``nbnd``\ (:math:`K^{el}`) 

Wave-number grid
----------------

- The :math:`{\bf k}` grid for the electronic state in :ref:`scf`, :ref:`ph`

   It is specified in the input file of ``pw.x`` as follows:
  
   ::
      
      K_POINTS automatic 
      {nk1} {nk2} {nk3} 0 0 0


   The numerical cost for :ref:`scf` and :ref:`ph` is proportional to :math:`N_{\bf k}^{\rm smooth}`
   (the number of :math:`{\bf k}` points in this grid).
      
- The :math:`{\bf q}` grid for :ref:`phonon`, the :math:`{\bf k}` grid for :ref:`twin`

   ``nq1, nq2, nq3`` in the input of ``ph.x``, 
   arguments of :ref:`twingrid`, and
   ``nk1, nk2, nk3`` in the input of :ref:`elph` must be the same.

   The :math:`N_{\bf q}` (the number of :math:`{\bf q}` in this grid) dependence of
   each program becomes as follows:

   - The numerical cost of ``pw.x`` in :ref:`twin` is proportional to :math:`N_{\bf q}`.
      
   - The numerical cost for all :math:`{\bf q}` in :ref:`ph` is proportional to :math:`N_{\bf q}`.

   - The numerical cost for all :math:`{\bf q}` in :ref:`elph` is proportional to :math:`N_{\bf q}^2`.

   - The numerical cost for all :math:`{\bf q}` in :ref:`sctk.x <sctk>` is proportional to :math:`N_{\bf q}^2`.

- The :math:`{\bf k}` grid in :ref:`dense` :ref:`[1] <ref>`

   In this calculation, the :math:`{\bf k}` grid should be as dense as that for the calculation
   of the density of states.
   The :math:`N_{\bf k}^{\rm dense}` (the number of :math:`{\bf k}` in this grid) dependence of
   each program becomes as follows:

   - Numerical costs for :ref:`scdft` and :ref:`sctk.x <sctk>` are not so affected by :math:`N_{\bf k}^{\rm dense}`.

   - The numerical cost of :ref:`deltaf` is proportional to :math:`N_{\bf k}^{\rm dense}`.

The relation of these :math:`{\bf k}` grid becomes as follows:

:math:`N_{\bf q} \leq N_{\bf k}^{\rm smooth} \leq N_{\bf k}^{\rm dense}`
