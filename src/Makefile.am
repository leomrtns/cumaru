# "-I$(srcdir)/.." means that I must "include <WFA/gap_affine/>" etc (i.e. the WFA dir is needed)
# "-I$(srcdir)/../kalign"  means that "include <kalign.h>" works
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../WFA -I$(srcdir)/../kalign -I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  ../WFA/libwfa.la ../kalign/libkalign.la ../biomcmc-lib/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = cumaru 
cumaru_SOURCES = main.c kseq.h
cumaru_LDADD = $(LDADD)
