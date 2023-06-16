#!/bin/sh
#SBATCH -N 1
#SBATCH -n 24

MPIRUN_OPTION="-bind-to board"
np_pw=24
np_kel=4
np_scdft=1
nt_pw=1
nt_kel=6
nt_scdft=24

export OMP_NUM_THREADS=${nt_pw}
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/pw.x -nk ${np_pw} -in scf.in > scf.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/ph.x -nk ${np_pw} -in ph.in > ph.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/ph.x -nk ${np_pw} -in elph.in > elph.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/ph.x -nk ${np_pw} -in epmat.in > epmat.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/q2r.x -in q2r.in > q2r.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/matdyn.x -in disp.in > disp.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/matdyn.x -in phdos.in > phdos.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/pw.x -nk ${np_pw} -in band.in > band.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/bands.x -nk ${np_pw} -in bands.in > bands.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/pw.x -nk ${np_pw} -in nscf.in > nscf.out 2>&1
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/projwfc.x -nk ${np_pw} -in proj.in > proj.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/fermi_velocity.x -in nscf.in > vfermi.out 2>&1
mpirun ${MPIRUN_OPTION} -np 1 ../../../bin/fermi_proj.x -in proj.in > pfermi.out 2>&1
cp twin0.in twin.in
../../../bin/twingrid.x 4 4 4 >> twin.in
mpirun ${MPIRUN_OPTION} -np ${np_pw} ../../../bin/pw.x -nk ${np_pw} -in twin.in > twin.out 2>&1
sed -e "/calculation/c calculation = \"kel\"" sctk0.in > kel.in
export OMP_NUM_THREADS=${nt_kel}
mpirun ${MPIRUN_OPTION} -np ${np_kel} ../../../bin/sctk.x -nk ${np_kel} -in kel.in > kel.out 2>&1
export OMP_NUM_THREADS=${nt_scdft}
sed -e "/calculation/c calculation = \"lambda_mu_k\"" sctk0.in > lambda_mu_k.in
mpirun ${MPIRUN_OPTION} -np ${np_scdft} ../../../bin/sctk.x -nk ${np_scdft} -in lambda_mu_k.in > lambda_mu_k.out 2>&1
sed -e "/calculation/c calculation = \"scdft_tc\"" sctk0.in > tc.in
mpirun ${MPIRUN_OPTION} -np ${np_scdft} ../../../bin/sctk.x -nk ${np_scdft} -in tc.in > tc.out 2>&1
sed -e "/calculation/c calculation = \"deltaf\"" sctk0.in > deltaf.in
mpirun ${MPIRUN_OPTION} -np ${np_scdft} ../../../bin/sctk.x -nk ${np_scdft} -in deltaf.in > deltaf.out 2>&1
sed -e "/calculation/c calculation = \"qpdos\"" sctk0.in > qpdos.in
mpirun ${MPIRUN_OPTION} -np ${np_scdft} ../../../bin/sctk.x -nk ${np_scdft} -in qpdos.in > qpdos.out 2>&1
