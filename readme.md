# Installation

To use main branch

``` bash
$ git clone https://gitlab.com/QEF/q-e.git -b master
$ cd q-e
$ git clone https://github.com/mitsuaki1987/sctk.git -b main SCTK
$ patch -p1 < SCTK/patch.diff
```

To try the develop branch 

``` bash
$ git clone https://gitlab.com/QEF/q-e.git -b develop
$ cd q-e
$ git clone https://github.com/mitsuaki1987/sctk.git -b develop SCTK
$ patch -p1 < SCTK/patch.diff
```

Configure the environment with the script `configure`
as the same as the original Quantum ESPRESSO.
               
``` bash
$ ./configure --enable-openmp
$ make pw ph pp sctk
$ make
```

The executable file is `sctk.x`.
