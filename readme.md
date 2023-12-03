# Installation

To use main branch

``` bash
$ git clone https://gitlab.com/QEF/q-e.git
$ cd q-e
$ git checkout 96cdd5ac6af9c060be392a95f14dbcbca5c1a890
$ git clone https://github.com/mitsuaki1987/sctk.git -b main
$ patch -p1 < sctk/patch.diff
```

To try the develop branch 

``` bash
$ git clone https://gitlab.com/QEF/q-e.git
$ cd q-e
$ git checkout 408b83ffe041751e9db61dc64feb9f34d48a9556
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
