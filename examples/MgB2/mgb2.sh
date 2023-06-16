#!/bin/sh
#SBATCH -N 1
#SBATCH -n 24
MPIRUN_OPTION="-bind-to board"
export OMP_NUM_THREADS=1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/pw.x -nk 24 -in scf.in > scf.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/ph.x -nk 24 -in ph.in > ph.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/ph.x -nk 24 -in elph.in > elph.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/ph.x -nk 24 -in epmat.in > epmat.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/q2r.x -in q2r.in > q2r.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/matdyn.x -in disp.in > disp.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/matdyn.x -in phdos.in > phdos.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/pw.x -nk 24 -in band.in > band.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/bands.x -nk 24 -in bands.in > bands.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/pw.x -nk 24 -in nscf.in > nscf.out 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/projwfc.x -nk 24 -in proj.in > proj.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/fermi_velocity.x -in nscf.in > vfermi.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/fermi_proj.x -in proj.in > pfermi.out 2>&1
cp twin0.in twin.in
../../../bin/twingrid.x 4 4 3 >> twin.in
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/pw.x -nk 24 -in twin.in > twin.out 2>&1
sed -e "/calculation/c calculation = \"kel\"" sctk0.in > kel.in 2>&1
mpirun ${MPIRUN_OPTION} -np 24 ../../../bin/sctk.x -nk 24 -in kel.in > kel.out 2>&1
export OMP_NUM_THREADS=24
sed -e "/calculation/c calculation = \"lambda_mu_k\"" sctk0.in > lambda_mu_k.in 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/sctk.x -nk 1 -in lambda_mu_k.in > lambda_mu_k.out 2>&1
sed -e "/calculation/c calculation = \"scdft_tc\"" sctk0.in > tc.in 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/sctk.x -nk 1 -in tc.in > tc.out 2>&1
sed -e "/calculation/c calculation = \"deltaf\"" sctk0.in > deltaf.in 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/sctk.x -nk 1 -in deltaf.in > deltaf.out 2>&1
sed -e "/calculation/c calculation = \"qpdos\"" sctk0.in > qpdos.in 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/sctk.x -nk 1 -in qpdos.in > qpdos.out 2>&1
