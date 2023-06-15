.. _QuantumESPRESSO: https://www.quantum-espresso.org/resources/users-manual

Installation
============

Prerequisite
------------

The prerequisite of SCTK is the same as that of  QuantumESPRESSO_.

Installation procedure
----------------------

#. Clone thee source code of SCTK as follows:

   .. code-block:: bash

      $ git clone git://git.osdn.net/gitroot/sctk/sctk.git

   This was forked from the original QuantumESPRESSO_.

#. Configure the environment with the script ``configure``
   as the same as the original QuantumESPRESSO_.
               
   .. code-block:: bash

       $ ./configure

#. Make

   .. code-block:: bash

       $ make sctk

   The executable file is ``sctk.x``.
