# Makefile for Superconducting Toolkit

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PP/src \
          $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG).

SCTKOBJS = \
nesting.o \
nesting_proj.o \
sctk_broyden.o \
sctk_clock.o \
sctk_cnt_dsp.o \
sctk_coulomb.o \
sctk_dmuxc.o \
sctk_eliashberg_sub.o \
sctk_gapeq_rhs.o \
sctk_gauss_legendre.o \
sctk_ini_delta.o \
sctk_invert.o \
sctk_io_delta.o \
sctk_kernel_weight.o \
sctk_qpdos.o \
sctk_read_file.o \
sctk_rotate_kernel.o \
sctk_spinfluc.o \
sctk_stdin.o \
sctk_tetra.o \
sctk_usonic.o \
sctk_val.o \
sctk_wfc.o \
sctk_write_hash.o \
sctk_z.o

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
PHOBJS = ../../PHonon/PH/libph.a
PPOBJS = ../../PP/src/libpp.a
QEMODS = $(BASEMODS)
TLDEPS= pw ph pp

all : tldeps libs-sctk sctk.x twingrid.x nesting_proj.x nesting.x sctk_eliashberg.x 

libs-sctk : libsctk.a

hash:
	sed -e "s/HASH/`git show |awk 'NR==1{print $$2}'`/g" sctk_write_hash > sctk_write_hash.f90

sctk.x : sctk_main.o libsctk.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS) $(PHOBJS) $(PPOBJS)
	$(LD) $(LDFLAGS) -o $@ sctk_main.o libsctk.a \
		 $(PPOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/$@ . )

sctk_eliashberg.x : sctk_eliashberg.o libsctk.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS) $(PHOBJS) $(PPOBJS)
	$(LD) $(LDFLAGS) -o $@ sctk_eliashberg.o libsctk.a \
		 $(PPOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/$@ . )

nesting_proj.x : nesting_proj.o $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS) $(PHOBJS) $(PPOBJS)
	$(LD) $(LDFLAGS) -o $@ nesting_proj.o \
		 $(PPOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/$@ . )

nesting.x : nesting.o $(MODULES) $(PPOBJS)
	$(LD) $(LDFLAGS) -o $@ nesting.o $(PPOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/$@ . )

twingrid.x : twingrid.sh
	cp twingrid.sh twingrid.x
	chmod u+x twingrid.x
	( cd ../../bin ; ln -fs ../SCTK/src/twingrid.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libsctk.a : hash $(SCTKOBJS)
	$(AR) $(ARFLAGS) $@ $(SCTKOBJS)
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/sctk.x

include make.depend
# DO NOT DELETE
