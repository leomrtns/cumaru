AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@ @ZLIB_LIBS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @ZLIB_LIBS@ ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

common_headers = kalign.h \
tldevel.h \
rng.h \
global.h \
parameters.h \
alignment_parameters.h \
bisectingKmeans.h \
sequence_distance.h \
euclidean_dist.h \
alignment.h \
misc.h \
bpm.h

common_src = run_kalign.c \
tldevel.c \
rng.c \
parameters.c \
alignment_parameters.c \
bisectingKmeans.c \
sequence_distance.c \
euclidean_dist.c \
alignment.c \
misc.c \
bpm.c

noinst_LTLIBRARIES = libkalign.la   ## change to noinst_LT... if should be linked statically (not installed globally)
libkalign_la_SOURCES = config.h $(common_headers) $(common_src) 

