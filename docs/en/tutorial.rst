.. _FermiSurfer: http://fermisurfer.osdn.jp/
.. _pw.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PW.html
.. _ph.x: file:///C:/Users/kawamuura/program/qe/qe-dev/PW/Doc/INPUT_PH.html
.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

Tutorial
========

The following tutorial should be done in ``SCTK/examples/Al/``.

.. _scf:

SCF calculation of the charge density
-------------------------------------

**Input file**: scf.in

**Program**: pw.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ export OMP_NUM_THREADS=1
   $ mpiexec -np 29 PATH/pw.x -nk 29 -in scf.in > scf.out

**Important parameters**

    calculation = "scf"
        Self-consistent calculation is performed with pw.x_.

.. note::

   We can perform first whether :ref:`phonon` or `coulomb` .

.. _phonon:

Calculations of phonon and electron-phonon interaction
------------------------------------------------------

.. _ph:

Calculation of phonon frequency and deformation potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input file**: ph.in

**Program**: ph.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ mpiexec -np 29 PATH/ph.x -nk 29 -in ph.in > ph.out

**Important parameters**

    fildvscf = 'dv'
        The file name for the deformation potential.
        It will be used in the next step.

    ldisp = .true.
        We compute phonons on the uniform :math:`{\bf q}` grid.

    lshift_q = .true.
        We shift the :math:`{\bf q}` grid for avoiding
        the singularity at :math:`\Gamma` point.

    nq1, nq2, nq3
        The :math:`q` grid.

.. _elph:

Calculation of electron-phonon interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input file**: elph.in

**Program**: ph.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ mpiexec -np 8 PATH/ph.x -nk 8 -in epmat.in > epmat.out

**Important parameter**

    electron_phonon = "scdft_input"
        compute the electron-phonon vertex from the deformation potential and
        the dynamical matrix that are already computed.

    elph_nbnd_min, elph_nbnd_max
        Since the electron-phonon interaction between the electronic states
        only in the vicinity of the Fermi surface affect the gap equation,
        we can reduce bands to decrease the numerical cost.
        The upper- and the lower limit of the bands that contain the Fermi level
        can be obtained by a program
        `fermi_velocity.x <https://www.quantum-espresso.org/Doc/pp_user_guide/>`_
        in QuantumESPRESSO_.

.. _coulomb:
   
Calculation of screened Coulomb interaction/Spin-fluctuation
------------------------------------------------------------

.. _dense:

Non-SCF calculation with a dense k grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input file**: nscf.in

**Program**: pw.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ mpiexec -np 32 PATH/pw.x -nk 32 -in nscf.in > nscf.out

**Important parameter**

   calculation = "nscf"
       Perform non self-consistent calculation.
   
   la2f = .true.
       Generate a file pwscf.a2Fsave which contains Kohn-Sham energy.

   nbnd
       For the calculating of the polarization function,
       we have to compute some empty states.
       We do not have to include so many empty state
       as the calculation of the insulator. 
       Typically, the number of empty states becomes
       the same as the number of occupied states.

.. _twin:

Calculation of wave functions for the screened Coulomb interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input file**: twin.in

**Program**: pw.x_ (QuantumESPRESSO_)

.. code-block:: bash

   $ bash PATH/twingrid.x 4 4 3 >> twin.in
   $ mpiexec -np 32 PATH/pw.x -nk 32 -in twin.in > twin.out

**Important parameter**

   calculation = "bands"
       Generated :math:`k` points by using
       :ref:`twingrid` and redirect it to the input file as above.
       This :math:`k` mesh must be the same as nq1, nq2, and nq3
       in the input of ph.x_ for :ref:`elph`.
       Then run pw.x_ with this input file.

Calculation of screened Coulomb interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input file**: sctk.in

**Program**: :ref:`sctk.x <sctk>`

.. code-block:: bash

   $ mpiexec -np 32 PATH/sctk.x -nk 32 -in sctk.in > kel.out

**Important parameters**

    :ref:`calculation = "kel" <kel>`
         Calculation of screened Couplmb / spin-fluctuation interaction.

    nq1, nq2,  nq3
         They must be the same as the :math:`{\bf k}` in the previous step.
   
.. _scdftscf:
   
SCDFT SCF calculation
---------------------

After all calculation in :ref:`phonon` and :ref:`coulomb` finished,
we can start SCDFT calculation.

**Input file**: sctk.in (Should be modified)

**Program**: :ref:`sctk.x <sctk>`

.. code-block:: bash

   $ export OMP_NUM_THREADS=32
   $ mpiexec -np 1 PATH/sctk.x < sctk.in > tc.out

**Important parameters**

    :ref:`calculation = "scdft_tc" <scdfttc>`
         By changing this part, we change the type of calculation.
         Here, we obtain :math:`T_c` by the bisection method.

Further analysis
----------------

By changing this part, we can perform the following analysis.

**Input file**: sctk.in (Should be modified)

**Program**: :ref:`sctk.x <sctk>`

.. code-block:: bash

   $ mpiexec -np 1 PATH/sctk.x < sctk.in

**Important parameters**

    :ref:`calculation = "lambda_mu_k" <lambdamuk>`
        
        Output data for FermiSurfer_ to plot :math:`\lambda_{n {\bf {\bf k}}}`.

    :ref:`calculation = "scdft" <scdft>`
        Perform SCDFT calculation at a temperature,
        then output a file delta.dat which is used by the following
        post-processes.

    temp
       Temperature (Kelvin).
       If this is set to 0 or negative values,
       sctk.x employs the special algorithm for the exact zero Kelvein.
       
    :ref:`calculation = "deltaf" <deltaf>`
        Output data for FermiSurfer_ to plot :math:`\Delta_{n {\bf k}}`.

    :ref:`calculation = "qpdos" <qpdos>`
         Quasi particle DOS. It requires long computational time.

.. note::

   There is another tutorial in SCTK/examples/MgB2/.
   Please note that in these tutorials, the number of
   :math:`{\bf k}` points and bands and the pseudopotentials are
   not sufficient for the production level.
