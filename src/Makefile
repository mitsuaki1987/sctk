# Makefile for Superconducting Toolkit

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PP/src \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../Modules \
          $(MOD_FLAG).

IFLAGS=

SCTKOBJS = \
    nesting_proj.o \
    sctk_broyden.o \
    sctk_cnt_dsp.o \
    sctk_clock.o \
    sctk_coulomb.o \
    sctk_dmuxc.o \
    sctk_gapeq_rhs.o \
    sctk_gauss_legendre.o \
    sctk_ini_delta.o \
    sctk_invert.o \
    sctk_io_delta.o \
    sctk_kernel_weight.o \
    sctk_main.o \
    sctk_qpdos.o \
    sctk_read_file.o \
    sctk_rotate_kernel.o \
    sctk_stdin.o \
    sctk_spinfluc.o \
    sctk_tetra.o \
    sctk_usonic.o \
    sctk_val.o \
    sctk_wfc.o \
    sctk_write_hash.o \
    sctk_z.o

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
PHOBJS = ../../PHonon/PH/libph.a
PPOBJS = ../../PP/src/libpp.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../upflib/libupf.a ../../XClib/xc_lib.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a

TLDEPS= pw ph pp

all : tldeps libs-sctk sctk.x twingrid.x nesting_proj.x

libs-sctk : libsctk.a

hash:
	sed -e "s/HASH/`git show | awk 'NR==1{print $$2}'`/g" sctk_write_hash > sctk_write_hash.f90

sctk_z.o : sctk_z.f90
	$(F90) $(FFLAGS) $(MODFLAGS) -c sctk_z.f90 -o sctk_z.o

sctk_write_hash.o : sctk_write_hash.f90
	$(F90) $(FFLAGS) $(MODFLAGS) -c sctk_write_hash.f90 -o sctk_write_hash.o

sctk.x : sctk_main.o libsctk.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS) $(PHOBJS) $(PPOBJS)
	$(LD) $(LDFLAGS) -o $@ sctk_main.o libsctk.a \
		 $(PPOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/sctk.x . )

nesting_proj.x : nesting_proj.o $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS) $(PHOBJS) $(PPOBJS)
	$(LD) $(LDFLAGS) -o $@ nesting_proj.o \
		 $(PPOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/nesting_proj.x . )

twingrid.x : twingrid.sh
	cp twingrid.sh twingrid.x
	chmod u+x twingrid.x
	( cd ../../bin ; ln -fs ../SCTK/src/twingrid.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libsctk.a : hash sctk_write_hash.o sctk_z.o $(SCTKOBJS)
	$(AR) $(ARFLAGS) $@ $(SCTKOBJS) sctk_z.o
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/sctk.x

include make.depend
# DO NOT DELETE

