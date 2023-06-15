# Makefile for Superconducting Toolkit

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PP/src \
          $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG).
IFLAGS=

SCTKOBJS = \
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
sctk_z.o

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
PHOBJS = ../../PHonon/PH/libph.a
PPOBJS = ../../PP/src/fermisurfer_common.o
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a

TLDEPS= pw ph pp

all : tldeps libs-sctk sctk.x twingrid.x

libs-sctk : libsctk.a

sctk.x : sctk_main.o libsctk.a $(PWOBJS) $(LRMODS) $(PHOBJS) $(PPOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ sctk_main.o libsctk.a \
		 $(PPOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../SCTK/src/sctk.x . )

twingrid.x : twingrid.sh
	cp twingrid.sh twingrid.x
	chmod u+x twingrid.x
	( cd ../../bin ; ln -fs ../SCTK/src/twingrid.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libsctk.a : $(SCTKOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/sctk.x

include make.depend
# DO NOT DELETE
