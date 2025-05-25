# Installation

To use main branch

``` bash
$ git clone https://gitlab.com/QEF/q-e.git
$ cd q-e
$ git checkout qe-7.4.1
$ git clone https://github.com/mitsuaki1987/sctk.git -b main
$ patch -p1 < sctk/patch.diff
```

To try the develop branch 

``` bash
$ git clone https://gitlab.com/QEF/q-e.git
$ cd q-e
$ git checkout qe-7.4.1
$ git clone https://github.com/mitsuaki1987/sctk.git -b develop
$ patch -p1 < sctk/patch.diff
```

Configure the environment with the script `configure`
as the same as the original Quantum ESPRESSO.
               
``` bash
$ ./configure --enable-openmp
$ make pw ph pp
$ cd sctk
$ make
```

The executable file is `sctk.x`.
