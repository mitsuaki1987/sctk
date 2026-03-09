.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

Installation
============

Prerequisite
------------

The prerequisite of SCTK is the same as that of  QuantumESPRESSO_.

Installation procedure
----------------------

#. Clone thee source code of QuantumESPRESSO_  and SCTK as follows:

   .. code-block:: bash

      $ git clone https://gitlab.com/QEF/q-e.git
      $ cd q-e
      $ git checkout qe-7.4.1
      $ git clone https://github.com/mitsuaki1987/sctk.git -b develop SCTK
      $ patch -p1 < SCTK/patch.diff

   To try the developping version, the corresponding original QE hash can be seen
   at https://github.com/mitsuaki1987/sctk/blob/develop/readme.md  

   This was forked from the original QuantumESPRESSO_.

#. Configure the environment with the script ``configure``
   as the same as the original QuantumESPRESSO_.
               
   .. code-block:: bash

       $ ./configure

#. Make

   .. code-block:: bash

       $ make pw ph pp sctk

   The executable file is ``sctk.x``.
   Sometimes, we need to ``make`` again because of incomplete compilation of SCTK.

   .. code-block:: bash

       $ make sctk
