This subdirectory contains a contributed example by Alpin N. Tatan for the elemental Pb with production level parameters extracted from M. Kawamura, Y. Hizume, and T. Ozaki, Phys. Rev. B 101, 134511 (2020) made by modifying the original tutorial files on Al and $\mathrm{MgB}_2$.

Running these input files would require a supercomputing facility. 

The required nodes and parallelization conditions can be read from the output files.

The output values match closely the published value in the abovementioned paper: 

($T_c$ = 4.4 K and 3.7 K for SCDFT calculations without and with spin fluctuations, respectively.)

At this moment, the example files do not include the spin-orbit interaction.

The example contains three folders. 

Beside the folders for SCDFT calculations with/without spin fluctuations (SF), a folder containing a set of DFPT calculations leading to the $\lambda,\omega_\mathrm{ln}$ in the output file "lambda" of the matdyn.x program is provided. These parameters are used in the McMillan formula of superconducting $T_c$. 

The output values for these parameters are also quite close with published values in the abovementioned paper.

The order for running the files:

SCDFT:

1. scf.in       (pw.x)
2. ph.in        (ph.x)
3. epmat.in     (ph.x)
4. nscf.in      (pw.x)
5. twin.in      (pw.x)
6. kel.in       (sctk.x)
7. tc.in        (sctk.x)

McMillan:

1. scf.in
2. ph.in                     ---> steps 1-2 is identical with SCDFT.
3. elph.in       (ph.x)
4. epmat.in
5. q2r.in        (q2r.x)
6. disp.in       (matdyn.x)
7. phdos.in      (matdyn.x)  ---> this produces the file 'lambda' containing the McMillan parameters.

In addition, the vfermi.in can be run optionally after step 1, using scf.in as input to fermi_velocity.x to find the bands that contain the Fermi surface for the input of epmat.in
