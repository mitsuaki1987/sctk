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
      $ git checkout 96cdd5ac6af9c060be392a95f14dbcbca5c1a890
      $ git clone https://github.com/mitsuaki1987/sctk.git -b sctk1.2.1-qe6.7
      $ patch -p1 < sctk/patch.diff

   To try the developping version, the corresponding original QE hash can be seen
   at https://github.com/mitsuaki1987/sctk/blob/develop/readme.md  

   .. code-block:: bash

      $ git clone https://gitlab.com/QEF/q-e.git
      $ cd q-e
      $ git checkout THE_HASH_CAN_BE_SEEN_AT https://github.com/mitsuaki1987/sctk/blob/develop/readme.md
      $ git clone https://github.com/mitsuaki1987/sctk.git -b develop
      $ patch -p1 < sctk/patch.diff

   This was forked from the original QuantumESPRESSO_.

#. Configure the environment with the script ``configure``
   as the same as the original QuantumESPRESSO_.
               
   .. code-block:: bash

       $ ./configure

#. Make

   .. code-block:: bash

       $ make pw ph pp
       $ cd sctk
       $ make

   The executable file is ``sctk.x``.
